# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Boxplot on sum of arthritis score
# Data: Arthritis scores on CIA mice

# Load library ------------------------------------------------------------

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_blank(),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Load data ---------------------------------------------------------------

# Sum of arthritis score
score <- read.csv("./data/metadata/CIA_Ctrl_Score_Sum_Sacrifice_20210118.csv")
View(score)

# Timepoint match table
matchTb <- read.csv("./data/metadata/matchTable_absolutDay_Timepoint.csv")
View(matchTb)

# Tidy data ---------------------------------------------------------------

# Merge data
df <- matchTb %>%
  left_join(score, by = "absolutDay")
View(df)

# Define orders
df$Treatment <- factor(df$Treatment, levels = c("Ctrl", "CIA"))

df$Timepoint <- factor(df$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Boxplot with dots -------------------------------------------------------

basic_plot <- df %>%
  ggplot(aes(x = Timepoint, y = Score, fill = Treatment)) +
  geom_boxplot(
    position = position_dodge(0.9), outlier.shape = NA,
    alpha = .2, color = "grey55"
  ) +
  geom_point(
    size = .8, color = "grey30",
    position = position_jitterdodge(dodge.width = 0.9, seed = 1234)
  ) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 13)) +
  labs(y = "Arthritis score (Sum of 4 paws)", x = NULL) +
  scale_fill_brewer(palette = "Dark2")

# Add significance to plot
p <- basic_plot +
  annotate(geom = "text", x = 4, y = 12.5, label = "*", size = 3) +
  annotate(geom = "text", x = 5, y = 12.5, label = "**", size = 3) +
  annotate(geom = "text", x = 6, y = 12.5, label = "***", size = 3) +
  annotate(geom = "text", x = 7, y = 12.5, label = "***", size = 3) +
  annotate(geom = "text", x = 8, y = 12.5, label = "***", size = 3) +
  geom_segment(aes(x = 3.7, y = 12, xend = 4.3, yend = 12)) +
  geom_segment(aes(x = 4.7, y = 12, xend = 5.3, yend = 12)) +
  geom_segment(aes(x = 5.7, y = 12, xend = 6.3, yend = 12)) +
  geom_segment(aes(x = 6.7, y = 12, xend = 7.3, yend = 12)) +
  geom_segment(aes(x = 7.7, y = 12, xend = 8.3, yend = 12))

ggsave("./plot/phenotype/boxPlot_scoreSum.tif", plot = p, device = "tiff",
       units = "mm", width = 100, height = 75, dpi = 300)

# t-test ------------------------------------------------------------------

# Create a variable
pVal <- NULL

for (i in unique(df$Timepoint)) {
  
  df_subset <- df %>%
    filter(Timepoint == i)
  
  # Compute t-test
  res <- t.test(Score ~ Treatment, data = df_subset)
  
  # Print p value
  print(paste0(i, ": ", res$p.value))
  
  # Put into a data frame
  pVal <- rbind(pVal, data.frame(i, res$p.value))
  
}

View(pVal)

# Add stars according to p value
pVal <- pVal %>%
  mutate(Star = case_when(
    res.p.value > 0.01 & res.p.value < 0.05 ~ "*",
    res.p.value > 0.001 & res.p.value < 0.01 ~ "**",
    res.p.value < 0.001 ~ "***",
    TRUE ~ "Not Sig."
  ))

pVal

# Save p values
fwrite(pVal, "./result/Significance_Test/boxPlot_scoreSum_tTest.csv")
       