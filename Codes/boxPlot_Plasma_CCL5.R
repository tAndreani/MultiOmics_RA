# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Boxplot of Luminex results on CCL5
# Data: Luminex results, adjusted for dilution factor

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Import data -------------------------------------------------------------

# Luminex data
lumData <- fread("./data/luminex/Luminex_Score_Merge_All_20210118.csv", quote = "")
View(lumData)

# Batch data
batchData <- read.csv("./data/metadata/paw_sample_merged_all_batches.csv")
View(batchData)

# Tidy data ---------------------------------------------------------------

# Merge with timepoint information
dfMerge <- batchData %>%
  select(absolutDay, Timepoint, Treatment, mouseID) %>%
  distinct() %>%
  right_join(lumData, by = c("mouseID", "Treatment", "absolutDay")) %>%
  filter(absolutDay != -7)

View(dfMerge)

unique(dfMerge$absolutDay)
unique(dfMerge$Analyte)

# Rank levels
dfMerge$Treatment <- factor(dfMerge$Treatment, levels = c("Ctrl", "CIA"))
dfMerge$Timepoint <- factor(dfMerge$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Boxplot -----------------------------------------------------------------

basic_plot <- dfMerge %>%
  filter(Analyte == "CCL5/RANTES") %>%
  ggplot(aes(x = Timepoint, y = pg_ml, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.9), outlier.shape = NA, 
               alpha = .2, color = "grey45") +
  geom_point(aes(color = Score), size = .8, alpha = .8,
             position = position_jitterdodge(dodge.width = 0.9, seed = 1234)) +
  labs(y = "Plasma CCL5 (pg/ml)", x = NULL, color = "Score") +
  scale_y_continuous(breaks = seq(0, 300, 50), limits = c(0, 310)) +
  scale_colour_viridis_b(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_fill_brewer(palette = "Dark2")

# Add p values
p <- basic_plot +
  annotate(geom = "text", x = 1, y = 290, label = "***", size = 3) +
  geom_segment(aes(x = 0.7, y = 280, xend = 1.3, yend = 280)) 

ggsave("./plot/phenotype/boxPlot_ccl5.tif", plot = p, device = "tiff",
       units = "mm", width = 100, height = 60, dpi = 300)

# t-test ------------------------------------------------------------------

# Subset for analyte name
df <- dfMerge %>%
  filter(Analyte == "CCL5/RANTES")

# Create a variable
pVal <- NULL

for (i in unique(df$Timepoint)) {
  df_subset <- df %>%
    filter(Timepoint == i)
  
  # Compute t-test
  res <- t.test(pg_ml ~ Treatment, data = df_subset)
  
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
fwrite(pVal, "./result/Significance_Test/boxPlot_ccl5_tTest.csv")
