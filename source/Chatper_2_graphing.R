library(tidyverse)
library(latex2exp)

GOES <- read_csv("D:/main/Datasets/Solar_REU/GOES_data/goes_clean.csv")
GOES$cycle <- as.factor(GOES$cycle)
GOES <- GOES %>% filter(is.na(Gflrtotalenergy) == FALSE | Gflrtotalenergy > 0)
GOES$goessat <- factor(as.factor(GOES$goessat), levels = c("GO7", "GO8", "GO9", "G10", "G11", "G12", "G13", "G14", "G15", "G16"))

satellite_data <- GOES %>% select(goessat, Gpeak)
ggplot(data = satellite_data, aes(x = Gpeak, y = goessat, color = goessat)) + 
  geom_point() +
  theme_bw() +
  labs(x = "Year", y = "Active Satellite", color = "GOES Series")

ggplot(data = GOES, aes(x = log(Gflrtotalenergy, 10))) + 
  geom_histogram(bins = 1000) + 
  theme_bw() + 
  labs(y = "Count", x = TeX("$\\log_{10}$(Total Energy)"))
