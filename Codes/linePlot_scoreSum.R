# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Line plot of arthritis score development
# Data: Arthritis scores on CIA mice

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set plot parameter ------------------------------------------------------

theme_set(theme_classic())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_blank(),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm"),
  legend.position = "top"
)

# Load data ---------------------------------------------------------------

score <- fread("./data/metadata/CIA_Score_Sum.csv")
View(score)

# Tidy data ---------------------------------------------------------------

# Change to long format
score_long <- score %>%
  pivot_longer(
    cols = 4:147,
    names_to = "Mouse",
    values_to = "scoreSum",
    values_drop_na = TRUE
  )

View(score_long)

# Line plot: Calculate statistics -----------------------------------------

score_stat <- score_long %>%
  group_by(absolutDay) %>%
  summarise(
    avg = mean(scoreSum, na.rm = TRUE),
    sd = sd(scoreSum, na.rm = TRUE),
    count = n(),
    se = sd(scoreSum, na.rm = TRUE) / sqrt(n())
  )

View(score_stat)

# Define disease stages ---------------------------------------------------

# Define breaks for immunization
rects <- data.frame(
  xstart = c(0, 21), xend = c(21, 70),
  col = as.factor(c(
    "After 1st immunization", "After 2nd immunization"
  ))
)

rects

# Define arrows of sacrifice dates ----------------------------------------

d <- data.frame(
  x = c(2, 14, 23, 29, 42, 49, 56, 69), y = c(1.5, 2, 3, 4, 5, 6, 7, 8),
  vx = rep(0, 8), vy = rep(-1, 8)
)
d

# Line plot with immunization stages --------------------------------------

# Line plot with background colors
basic_plot <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.4) +
  geom_segment(
    data = d, mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy),
    arrow = arrow(length = unit(0.1, "cm")), size = 1, color = "gray50"
  ) +
  geom_line(data = score_stat, aes(x = absolutDay, y = avg), size = 1, color = "grey50") +
  geom_point(data = score_stat, aes(x = absolutDay, y = avg), size = 1.5, color = "grey4") +
  geom_errorbar(data = score_stat, aes(x = absolutDay, y = avg, ymin = avg - se, ymax = avg + se), size = .8, width = 0, color = "grey4", alpha = .4) +
  labs(x = "Day", y = "Arthritis score (Sum of 4 paws)", caption = "Mean ± SEM; Arrow: Sacrifice day", fill = NULL) +
  ylim(0, 8.5) +
  scale_x_continuous(breaks = seq(0, 70, 10)) +
  scale_fill_brewer(palette = "Pastel1")

# Add annotations for time points
p <- basic_plot +
  annotate(geom = "text", x = 2, y = 2, label = "2 days", size = 2.5) +
  annotate(geom = "text", x = 14, y = 2.5, label = "2 weeks", size = 2.5) +
  annotate(geom = "text", x = 23, y = 3.5, label = "3 weeks", size = 2.5) +
  annotate(geom = "text", x = 29, y = 4.5, label = "4 weeks", size = 2.5) +
  annotate(geom = "text", x = 42, y = 5.5, label = "6 weeks", size = 2.5) +
  annotate(geom = "text", x = 49, y = 6.5, label = "7 weeks", size = 2.5) +
  annotate(geom = "text", x = 56, y = 7.5, label = "8 weeks", size = 2.5) +
  annotate(geom = "text", x = 65, y = 8.5, label = "10 weeks", size = 2.5)

ggsave("./plot/phenotype/linePlot_scoreSum_annotation.tif",
  plot = p, device = "tiff",
  units = "mm", width = 88, height = 88, dpi = 300
)
