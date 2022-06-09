library(tidyverse)
library(wesanderson)
library(car)
library(VGAM)
library(qqplotr)
library(lubridate)

GOES <- read_csv("D:/main/Datasets/Solar_REU/GOES_data/goes_clean.csv")
GOES$cycle <- as.factor(GOES$cycle)

# butterfly plot
ggplot(data = GOES %>% filter(Glat != 0)  %>% filter(Glat > -55), aes(x = Gstart, y = Glat, color = cycle)) + 
  geom_point(size = 1.7) +
  theme_bw() + 
  labs(x = "Year", y = "Latitude") +
  scale_color_brewer(palette = "Blues")

# preparing histograms for flare frequencies
dates <- c(as.Date("1996-07-30 00:00:01 UTC"), 
           as.Date("2008-12-01 00:00:01 UTC"), 
           as.Date("2001-11-15 00:00:01 UTC"), 
           as.Date("2014-04-15 00:00:01 UTC"))
extrema <- c("Minimum", "Minimum", "Maximum", "Maximum")
solar_extrema <- data.frame(dates, extrema)

# aggregated flare frequency histogram
ggplot(data = GOES, aes(x=Gtstart, fill = cycle)) + 
  geom_histogram(center = 1996.5, binwidth = 1, color = "black") +
  geom_vline(data = solar_extrema, mapping = aes(xintercept = decimal_date(dates), linetype = extrema), size = .85) +
  labs(x = "Year", y = "Count") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues")

# aggregated flare frequency bar plot
ggplot(data = GOES, aes(x=year, fill = as.factor(cycle))) + 
  geom_bar(color = "black") +
  geom_vline(data = solar_extrema, mapping = aes(xintercept = decimal_date(dates), linetype = extrema), size = .85) +
  labs(x = "Year", y = "Count") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues")


# histogram of flare counts by class
ggplot(data = GOES %>% mutate(class= recode(class, "c('A', 'B')='A-B'")), aes(x=Gtstart, fill = cycle)) +
  geom_histogram(center = 1996.5, binwidth = 1, color = "black") +
  geom_vline(data = solar_extrema, mapping = aes(xintercept = decimal_date(dates), linetype = extrema), size = .85) +
  facet_wrap(~ class) +
  scale_y_log10() +
  labs(y = "Count", x = "Year") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues")

# aggregated flare frequency bar plot by class
ggplot(data = GOES %>% mutate(class= recode(class, "c('A', 'B')='A-B'")), aes(x=year, fill = as.factor(cycle))) + 
  geom_bar(color = "black") +
  geom_vline(data = solar_extrema, mapping = aes(xintercept = decimal_date(dates), linetype = extrema), size = .85) +
  facet_wrap(~class) +
  scale_y_log10() +
  labs(x = "Year", y = "Count") +
  theme_bw() +
  scale_fill_brewer(palette = "Blues")

# monthly flare counts
GOES_counts <- GOES %>% mutate(yearmonth = as.Date(as.yearmon(format(Gstart, format = "%Y-%m")))) %>% count(yearmonth, cycle)

ggplot(data = GOES_counts, aes(x = yearmonth, y = n, color = cycle)) + 
  geom_point(size = 2) +
  theme_bw() + 
  labs(x = "Year", y = "Count") +
  scale_color_brewer(palette = "Blues")



GOES_2013 <- GOES %>% filter(Garreg > 0) %>% filter(year == 2013) %>% 
  select(Garreg, year, cycle, Gpeak, Gstart, Gstop, Glat) %>% 
  group_by(Garreg) %>% 
  mutate(waiting_time = as.numeric(Gstart-lag(Gstart))) %>% 
  drop_na() %>% 
  summarize(start = min(Gstart), end = max(Gstop), y = mean(Glat[Glat != 0]), y1 = ifelse(is.na(y), 0, y), lambda = 1/mean(waiting_time))

ggplot(data = GOES_2013) + 
  geom_segment(aes(x = start, xend = end, y = y1, yend = y1, color = lambda), size = 1.1) +
  theme_bw() +
  labs(x = "Date", y = "Latitude") + scale_color_distiller(palette = "Blues")






# Removed first flare in active region for PP plot
AR_2403 <- GOES %>% filter(Garreg == AR_nums[2]) %>% select(Gstart, Gpeak, Gstop, Garreg, class, year) %>% mutate(waiting_time = as.numeric(Gstart-lag(Gstart))) %>% filter(is.na(waiting_time) == FALSE)

AR_num <- AR_nums[2]
AR_year <- AR_2403$year[1]
flare_count <- nrow(AR_2403) + 1
AR_length_days <- round(as.numeric(AR_2403$Gstart[length(AR_2403$Gstart)] - AR_2403$Gstart[1]),2)
sample_mean <- mean(AR_2403$waiting_time)
sample_sd <- sd(AR_2403$waiting_time)
ratio <- sample_sd/sample_mean
lambda <- 1/mean(as.numeric(AR_2403$waiting_time))
KS_pval <- ks.test(AR_2403$waiting_time, "pexp", rate = lambda)$p.value
results <- data.frame(AR_num, AR_year, flare_count, AR_length_days, sample_mean, sample_sd, ratio, lambda, KS_pval)
colnames(results) <- c("AR Number", "AR Year", "Flare Count", "AR Length (Days)", "Sample Mean", "Sample SD", "Ratio", "Lambda-hat", "KS p-value")

ggplot(AR_2403, aes(sample = waiting_time)) + 
  stat_pp_point(distribution = "exp",
                color = "lightblue") +
  stat_pp_line(distribution = "exp",
               color = "red") +
  theme_bw() +
  labs(x = "Theoretical Cumulative Probability" , 
       y = "Empirical Cumulative Probability",
       title = paste("Exponential P-P Plot for Observed Waiting Times in AR", AR_num), 
       subtitle = paste("Each of the", 
                        nrow(AR_2403), 
                        "points represent the waiting time in minutes between consecutive flares \nin AR", 
                        AR_num, "in", AR_year," which lasted", AR_length_days, "days."))









GOES_ar_count <- GOES_select %>% mutate(yearmonth = as.Date(as.yearmon(format(Gpeak, format = "%Y-%m")))) %>% count(yearmonth, cycle)

# Does not include Active Region 0
ggplot(data = GOES_ar_count, aes(x = yearmonth, y = n, color = cycle)) + 
  geom_point(size = 2) +
  theme_bw() + 
  labs(x = "Year", y = "Count") +
  scale_color_brewer(palette = "Blues")

GOES_select <- GOES_AR %>% select(Garreg, year, cycle, Gpeak, Glat) %>% group_by(Garreg) %>% filter(row_number() == 1) %>% ungroup()

# removed all ARs with mean lattitude of 0
GOES_ar_lat <- GOES_select %>% 
  mutate(yearmonth = as.Date(as.yearmon(format(Gtart, format = "%Y-%m")))) %>%
  group_by(Garreg) %>% 
  mutate(mean_lat = mean(Glat)) %>% 
  filter(mean_lat != 0)

ggplot(GOES_ar_lat %>% filter(Glat > -55), aes(x = yearmonth, y = mean_lat, color = cycle)) + 
  geom_point(size = 2) +
  theme_bw() + 
  scale_color_brewer(palette = "Blues") +
  labs(title = "Average Active Region Latitudes", x = "Year", y = "Average Latitude", subtitle = "Filtered out Active Regions with Average Latitudes of 0")








####
# Evaluating Goodness of Fit
####

##
# Probability - Probability Plots
##

# random sample of sample size 100 from exponential(5)
set.seed(10)
exp_sample <- rexp(100, rate = 5)
sim_data <- data.frame(sample = exp_sample)
lambda <- 1/mean(exp_sample)

ggplot(sim_data, aes(sample = sample)) +
  stat_pp_band(distribution = "exp") +
  stat_pp_line(color = "firebrick3") +
  stat_pp_point(distribution = "exp",
                color = "deepskyblue4") +
  theme_bw() +
  labs(y = "Empirical Cumulative Probability", x = "Theoretical Cumulative Probability")

set.seed(10)
log_sample <- rlnorm(100, meanlog = 0, sdlog = 2)
sim_data <- data.frame(sample = log_sample)
lambda <- 1/mean(log_sample)

ggplot(sim_data, aes(sample = sample)) + 
  stat_pp_band(distribution = "exp") +
  stat_pp_line(color = "firebrick3") +
  stat_pp_point(distribution = "exp",
                color = "deepskyblue4") +
  theme_bw() +
  labs(y = "Empirical Cumulative Probability", x = "Theoretical Cumulative Probability")

