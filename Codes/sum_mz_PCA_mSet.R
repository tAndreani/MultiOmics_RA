# Introduction ------------------------------------------------------------

# Image: metabolomicsr
# Output: PCA plot
# Data: Saved mSet object of MetaboAnalyst

# Load library ------------------------------------------------------------

library(data.table)
library(tidyverse)
library(MetaboAnalystR)
library(RColorBrewer)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Load mSet object --------------------------------------------------------

mSet <- readRDS(file = "./data/DMD1659/sum_mz_mSet_MetaboAnalyst_norm.rds")

# Calculate PCA -----------------------------------------------------------

mSet <- PCA.Anal(mSet)

# Package-contained command to generate PCA plot --------------------------

mSet <- PlotPCA2DScore(mSet, "PCA_2D_", "png", 72, width = 12, pcx = 1, pcy = 2, reg = 0.95, show = 1, grey.scale = 0)

# Create data frame with PCA scores ---------------------------------------

# Extract PC scores
score <- as.data.frame(mSet$analSet$pca$x)
View(score)

# Add meta data to PCA scores
score$Treatment <- mSet$dataSet$facA
score$Time <- mSet$dataSet$facB

# Ranks
score$Time <- factor(score$Time, levels = c("2 weeks", "4 weeks", "6 weeks", "10 weeks"))
score$Treatment <- factor(score$Treatment, levels = c("Ctrl", "CIA"))

# Define label names
xlabel <- paste0("PC1 ", "(", round(100 * mSet$analSet$pca$variance[1], 1), "%)")
xlabel # PC1 (19.2%)
ylabel <- paste0("PC2 ", "(", round(100 * mSet$analSet$pca$variance[2], 1), "%)")
ylabel # PC2 (17.7%)

# PCA with timepoints -----------------------------------------------------

# Create data frame to extract convex hull points
pca_hull <-
  score %>%
  group_by(Treatment) %>%
  slice(chull(PC1, PC2))

# PCA plot with convex hulls
pca_plot <- ggplot(score, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = Time, shape = Treatment), size = 2) +
  labs(x = xlabel, y = ylabel) +
  geom_polygon(data = pca_hull, aes(fill = Treatment), alpha = 0.1, show.legend = TRUE) +
  scale_shape_manual(values = c(1, 16)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_brewer(palette = "Dark2")

ggsave(filename = "./plot/DMD1659/sum_mz/PCA_timepoint_manuscript.tif", pca_plot, 
       device = "tiff", units = "mm", width = 92, height = 70, dpi = 300)

# Add frontleft score to PCA data -----------------------------------------

# Individual paw score at sacrifice
Score_Ind <- read.csv("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/metadata/CIA_Score_Individual_Sacrifice.csv", row.names = 1)
View(Score_Ind)

# Add FrontLeft score to PCA data
score_fl <- score %>%
  rownames_to_column(var = "mouseID") %>%
  mutate(mouseID = as.numeric(mouseID)) %>%
  left_join(Score_Ind[Score_Ind$Location == "FrontLeft", ], by = c("mouseID", "Treatment")) %>%
  mutate(Score = as.factor(replace_na(Score, 0)))
View(score_fl)

score_fl$Treatment <- factor(score_fl$Treatment, levels = c("Ctrl", "CIA"))

# PCA with FrontLeft score ------------------------------------------------

# Create data frame to extract convex hull points
pca_fl_hull <-
  score_fl %>%
  group_by(Treatment) %>%
  slice(chull(PC1, PC2))

# PCA with FrontLeft score
pca_plot2 <- score_fl %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = Score, shape = Treatment), size = 2) +
  labs(x = xlabel, y = ylabel, color = "Paw score") +
  geom_polygon(data = pca_fl_hull, aes(fill = Treatment), alpha = 0.1, show.legend = TRUE) +
  scale_shape_manual(values = c(1, 16)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02"))

ggsave(filename = "./plot/DMD1659/sum_mz/PCA_score_manuscript.tif", pca_plot2,
       device = "tiff", units = "mm", width = 90, height = 70, dpi = 300)
