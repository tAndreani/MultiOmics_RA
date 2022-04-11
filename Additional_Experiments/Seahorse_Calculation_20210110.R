# Introduction ------------------------------------------------------------

# Project: CIA (09/2020 - 12/2020)
# Experiment: Seahorse assay
# This R script produces violin plot and scatter plot

# Load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
cat("\f") # Clear console

# Import data -------------------------------------------------------------

# Data from all timepoints, except for 01102020
df_Pre <- fread("./Seahorse/data/Seahorse_tidyFormat_11242020.csv", quote = "")
View(df_Pre)

# Data from last sacrifice 01102020
df_01102020 <- fread("./Seahorse/data/Seahorse_01102020.csv", quote = "")
View(df_01102020)

# Import metadata
metaData <- fread("./Phenotype/data/CIA_spleenWeight.csv") %>%
  select(absolutDay, Treatment, mouseID) %>%
  as.data.frame()
View(metaData)

# Tidy --------------------------------------------------------------------

# Change to long format
a <- df_Pre %>%
  pivot_longer(
    cols = 5:19,
    names_to = "Time",
    values_to = "Value",
    values_drop_na = TRUE
  ) %>%
  select(!c(Timepoint, Sacrifice_Date))

b <- df_01102020 %>%
  pivot_longer(
    cols = 3:23,
    names_to = "mouseID",
    values_to = "Value",
    values_drop_na = TRUE
  )

# Merge
ab <- rbind(a, b[, c(1, 3, 2, 4)])
View(ab)

metaData$mouseID <- as.character(metaData$mouseID)

# Add metadata
df <- ab %>%
  left_join(metaData, by = "mouseID")
View(df)
glimpse(df)

# Tidy data
df$Time <- as.numeric(df$Time)
df$Parameter <- factor(df$Parameter, levels = c("OCR (pmol/min)", "ECAR (mpH/min)"))
df$Treatment <- factor(df$Treatment, levels = c("Ctrl", "CIA"))
df$mouseID <- as.factor(df$mouseID)
df$absolutDay <- as.factor(df$absolutDay)

# Save data
fwrite(df, file = "./Seahorse/data/Seahorse_Combined_All_20210110.csv")

# OCR and ECAR in one plot ------------------------------------------------

# Overview except for timepoint -7
df %>%
  filter(absolutDay != -7) %>%
  ggplot(aes(x = Time, y = Value, group = mouseID)) +
  geom_path(aes(color = Treatment), size = 1) +
  geom_point(size = .5, alpha = .7) +
  facet_grid(rows = vars(Parameter), cols = vars(absolutDay), scales = "free_y") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") ## Save as pdf

# Overview after removal of abnormal measurements
df %>%
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  ggplot(aes(x = Time, y = Value, group = mouseID)) +
  geom_path(aes(color = Treatment), size = 1) +
  geom_point(size = .5, alpha = .7) +
  facet_grid(rows = vars(Parameter), cols = vars(absolutDay), scales = "free_y") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")

# Mean with sd after removal of abnormal measurements
df_summmary <- df %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  group_by(Parameter, Time, absolutDay, Treatment) %>%
  summarise(len = mean(Value), sd = sd(Value))

p_summary <- df_summmary %>%
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  ggplot(aes(x = Time, y = len, group = Treatment, color = Treatment)) +
  geom_line() +
  geom_pointrange(aes(ymin = len - sd, ymax = len + sd), fatten = .2, size = 1) +
  labs(y = "Mean with standard deviation") +
  facet_grid(rows = vars(Parameter), cols = vars(absolutDay), scales = "free_y") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2")

ggsave("Mean_sd_20210110.png", plot = p_summary, width = 30, height = 10, units = "cm")

# Example plot
p_example <- df_summmary %>%
  filter(absolutDay == 2 & Treatment == "CIA") %>%
  ggplot(aes(x = Time, y = len)) +
  geom_line(aes(color = Parameter), size = 1) +
  geom_point(size = 2) +
  labs(y = NULL) +
  facet_grid(rows = vars(Parameter), scales = "free_y") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")

ggsave("Demo_OCR_ECAR_20210110.png", plot = p_example, width = 8, height = 8, units = "cm")

# Calculation -------------------------------------------------------------

# Non-Mitochondiral Respiration
## The median of the 3 OCR data points after antimycin A and Rotenone addition
df_nMitoR <- df %>%
  filter(Parameter == "OCR (pmol/min)" & Time %in% c(81, 87.7, 94.3)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(nMitoR = median(Value))
View(df_nMitoR)

# Basel OCR
## The mean OCR of the last 3 baseline data points minus non-mitochondrial respiration
df_BR <- df %>%
  filter(Parameter == "OCR (pmol/min)" & Time %in% c(1.3, 7.9, 14.5)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(Avg_First3 = mean(Value)) %>%
  left_join(df_nMitoR, by = c("absolutDay", "Treatment", "mouseID")) %>%
  mutate(BR = Avg_First3 - nMitoR)
View(df_BR)

# Proton Leak
## The mean of the 3 OCR data points after oligomycin addition minus non-mitochondrial respiration
df_PL <- df %>%
  filter(Parameter == "OCR (pmol/min)" & Time %in% c(41.2, 47.8, 54.4)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(Avg_p789 = mean(Value)) %>%
  left_join(df_BR, by = c("absolutDay", "Treatment", "mouseID")) %>%
  mutate(PL = Avg_p789 - nMitoR)
View(df_PL)

# ATP Production
## OCR (Before oligomycin) - OCR (after oligomycin)
df_ATP <- df %>%
  filter(Parameter == "OCR (pmol/min)" & Time %in% c(21.3, 27.9, 34.5)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(Avg_p456 = mean(Value)) %>%
  left_join(df_PL, by = c("absolutDay", "Treatment", "mouseID")) %>%
  mutate(ATP = Avg_p456 - PL)
View(df_ATP)

# Maximal Respiration
## The highest of the 3 OCR values after the addition of FCCP
df_MR <- df %>%
  filter(Parameter == "OCR (pmol/min)" & Time %in% c(61.1, 67.7, 74.3)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(Max_FCCP = max(Value)) %>%
  left_join(df_ATP, by = c("absolutDay", "Treatment", "mouseID")) %>%
  mutate(MR = Max_FCCP - nMitoR)
View(df_MR)

# Spare Respiratory Capacity
## Maximal Respiration - Basel OCR
df_SRC <- df_MR %>%
  mutate(SRC = MR - BR)

# Non-Glycolytic Acidification
## Basal ECAR
df_nGA <- df %>%
  filter(Parameter == "ECAR (mpH/min)" & Time %in% c(1.3, 7.9, 14.5)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(nGA = median(Value)) %>%
  left_join(df_SRC, by = c("absolutDay", "Treatment", "mouseID"))

# Glycolysis
df_Glc <- df %>%
  filter(Parameter == "ECAR (mpH/min)" & Time %in% c(21.3, 27.9, 34.5)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(AvgE_p456 = max(Value)) %>%
  left_join(df_nGA, by = c("absolutDay", "Treatment", "mouseID")) %>%
  mutate(Glc = AvgE_p456 - nGA)
View(df_Glc)

# Extra OCR and ECAR calculation
df_AvgE1 <- df %>%
  filter(Parameter == "ECAR (mpH/min)" & Time %in% c(1.3, 7.9, 14.5)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(AvgE_First3 = mean(Value)) %>%
  left_join(df_Glc, by = c("absolutDay", "Treatment", "mouseID"))

df_AvgE2 <- df %>%
  filter(Parameter == "ECAR (mpH/min)" & Time %in% c(61.1, 67.7, 74.3)) %>%
  group_by(absolutDay, Treatment, mouseID) %>%
  summarise(MaxE_FCCP = max(Value)) %>%
  left_join(df_AvgE1, by = c("absolutDay", "Treatment", "mouseID"))
View(df_AvgE2)

# Save data
fwrite(df_AvgE2, file = "./Seahorse/data/Seahorse_Calculation_All_20210110.csv")

# Plot calculation --------------------------------------------------------

# Tidy data
df_long <- df_AvgE2 %>%
  select(absolutDay, Treatment, mouseID, BR, ATP, MR, SRC, Glc) %>%
  pivot_longer(
    cols = 4:8,
    names_to = "Parameter",
    values_to = "Value",
    values_drop_na = TRUE
  )
View(df_long)

# Rank parameters
df_long$Parameter <- factor(df_long$Parameter, levels = c("BR", "Glc", "ATP", "MR", "SRC"))

# Violin plot
labelNames <- c(
  "BR" = "Basal respiration (pmol/min)",
  "ATP" = "ATP-coupled respiration (pmol/min)",
  "MR" = "Maximal respiration (pmol/min)",
  "SRC" = "Spare respiratory capacity (pmol/min)",
  "Glc" = "Glycolysis (mpH/min)"
)

p_cal <- df_long %>%
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  filter(Parameter != "SRC") %>%
  ggplot(aes(x = absolutDay, y = Value, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = .5, scale = "width", draw_quantiles = 0.5, size = 1, color = "#999999") +
  geom_point(position = position_jitterdodge(dodge.width = 0.9, seed = 2019), size = 1.5) +
  labs(x = "Day after immunization", y = NULL) +
  facet_grid(rows = vars(Parameter), scales = "free_y", labeller = as_labeller(labelNames)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"))

ggsave("Calculation_20210110.png", plot = p_cal, width = 25, height = 27, units = "cm")

# Plot OCR:ECAR -----------------------------------------------------------

# Basal
p_basal <- df_AvgE2 %>% 
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  ggplot(aes(x = AvgE_First3, y = Avg_First3)) +
  geom_point(aes(shape = Treatment, color = absolutDay), size = 2) +
  labs(x = "Basal ECAR (mpH/min)", y = "Basal OCR (pmol/min)") +
  scale_shape_manual(values = c(19, 21)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_viridis_d() +
  facet_grid(rows = vars(absolutDay))

ggsave("Basal_OCR_ECAR_20210110.png", plot = p_basal, width = 8, height = 10, units = "cm")

# Maximal
p_maximal <- df_AvgE2 %>% 
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  ggplot(aes(x = MaxE_FCCP, y = Max_FCCP)) +
  geom_point(aes(shape = Treatment, color = absolutDay), size = 2) +
  labs(x = "Maximal ECAR (mpH/min)", y = "Maximal OCR (pmol/min)", color = "Day after immunization") +
  scale_shape_manual(values = c(19, 21)) +
  theme_bw() +
  scale_color_viridis_d() +
  facet_grid(rows = vars(absolutDay))

ggsave("Maximal_OCR_ECAR_20210110.png", plot = p_maximal, width = 12.75, height = 10, units = "cm")
