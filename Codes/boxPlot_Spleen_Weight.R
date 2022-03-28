# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Boxplot on spleen weight
# Data: CIA mouse metadata

# Load library ------------------------------------------------------------

library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(qvalue)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Load data ---------------------------------------------------------------

# Batch data
batchData <- read.csv("./data/metadata/paw_sample_merged_all_batches.csv")
View(batchData)
glimpse(batchData)

# Sum of arthritis score
score <- read.csv("./data/metadata/CIA_Ctrl_Score_Sum_Sacrifice_20210118.csv")
View(score)

# Tidy data ---------------------------------------------------------------

# Merge data
df <- batchData %>%
  left_join(score, by = c("mouseID", "Treatment", "absolutDay"))
View(df)

# Define orders
df$Treatment <- factor(df$Treatment, levels = c("Ctrl", "CIA"))

df$Timepoint <- factor(df$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Boxplot with dots -------------------------------------------------------

basic_plot <- df %>%
  ggplot(aes(x = Timepoint, y = Spleen_mg, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA, 
               alpha = .2, color = "grey45") +
  geom_point(aes(color = Score), size = .8, alpha = .8,
             position = position_jitterdodge(dodge.width = 0.9, seed = 1234)) +
  labs(y = "Spleen weight (mg)", x = NULL, color = "Score") +
  scale_y_continuous(breaks = seq(50, 300, 50), limits = c(40, 325)) +
  scale_colour_viridis_b(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_fill_brewer(palette = "Dark2")

p <- basic_plot +
  annotate(geom = "text", x = 2, y = 320, label = "***", size = 3) +
  annotate(geom = "text", x = 3, y = 320, label = "***", size = 3) +
  annotate(geom = "text", x = 4, y = 320, label = "**", size = 3) +
  annotate(geom = "text", x = 5, y = 320, label = "***", size = 3) +
  annotate(geom = "text", x = 7, y = 320, label = "**", size = 3) +
  annotate(geom = "text", x = 8, y = 320, label = "**", size = 3) +
  geom_segment(aes(x = 1.7, y = 310, xend = 2.3, yend = 310)) +
  geom_segment(aes(x = 2.7, y = 310, xend = 3.3, yend = 310)) +
  geom_segment(aes(x = 3.7, y = 310, xend = 4.3, yend = 310)) +
  geom_segment(aes(x = 4.7, y = 310, xend = 5.3, yend = 310)) +
  geom_segment(aes(x = 6.7, y = 310, xend = 7.3, yend = 310)) +
  geom_segment(aes(x = 7.7, y = 310, xend = 8.3, yend = 310)) 

ggsave("./plot/phenotype/boxPlot_spleenWeight.tif", plot = p, device = "tiff",
       units = "mm", width = 100, height = 60, dpi = 300)

# t-test ------------------------------------------------------------------

# Create a variable
pVal <- NULL

for (i in unique(df$Timepoint)) {
  
  df_subset <- df %>%
    filter(Timepoint == i)
  
  # Compute t-test
  res <- t.test(Spleen_mg ~ Treatment, data = df_subset)
  
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
fwrite(pVal, "./result/Significance_Test/boxPlot_spleenWeight_tTest.csv")
