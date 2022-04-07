# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Boxplot with dots on m/z 664.11668 NAD
# Data: tMSI raw data

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 8),
  plot.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.key.size = unit(4, "mm")
)

# Load data ---------------------------------------------------------------

# Load tidied-up annotation file
annoMol.metaID <- fread("./data/DMD1659/DMD1659_ind_mz_annotation_tidy.csv")
View(annoMol.metaID)

# Load tidied-up raw data
rawData <- fread("./data/DMD1659/DMD1659_ind_mz_v2_raw_tidy.csv")
View(rawData)

# Individual paw score at sacrifice
Score_Ind <- read.csv("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/metadata/CIA_Score_Individual_Sacrifice.csv", row.names = 1)
View(Score_Ind)

# Subset for FrontLeft score
Score_FrontLeft <- Score_Ind %>%
  filter(Location == "FrontLeft") %>%
  select(Score, mouseID, Treatment)
View(Score_FrontLeft)

# Tidy data ---------------------------------------------------------------

# Extract X664.1167
dataTable <- rawData %>%
  select(mouseID, Treatment, Time, X664.1167) %>%
  left_join(Score_FrontLeft, by = c("mouseID", "Treatment")) %>%
  mutate(Score = replace(Score, is.na(Score), 0))
View(dataTable)

# Define orders
dataTable$Treatment <- factor(dataTable$Treatment, levels = c("Ctrl", "CIA"))

dataTable$Time <- factor(dataTable$Time,
  levels = c("2 weeks", "4 weeks", "6 weeks", "10 weeks")
)

# Boxplot with dots -------------------------------------------------------

basic_plot <- dataTable %>%
  ggplot(aes(x = Time, y = X664.1167, fill = Treatment)) +
  geom_boxplot(
    position = position_dodge(0.9), outlier.shape = NA,
    alpha = .2, color = "grey45"
  ) +
  geom_point(aes(color = Score),
    size = .8, alpha = .8,
    position = position_jitterdodge(dodge.width = 0.9, seed = 1234)
  ) +
  labs(y = "Signal quantification (664.11668 m/z)", x = NULL, color = "Paw score") +
  scale_colour_viridis_b(limits = c(0, 4), breaks = c(1, 2, 3)) +
  scale_fill_brewer(palette = "Dark2") +
  ylim(200000, 800000)

# Add p values
p <- basic_plot +
  annotate(geom = "text", x = 1, y = 790000, label = "**", size = 3) +
  annotate(geom = "text", x = 2, y = 790000, label = "*", size = 3) +
  annotate(geom = "text", x = 3, y = 790000, label = "*", size = 3) +
  geom_segment(aes(x = 0.7, y = 770000, xend = 1.3, yend = 770000)) +
  geom_segment(aes(x = 1.7, y = 770000, xend = 2.3, yend = 770000)) +
  geom_segment(aes(x = 2.7, y = 770000, xend = 3.3, yend = 770000))

ggsave("./plot/DMD1659/ind_mz/boxPlot_NAD_Manuscript.tif",
  plot = p, device = "tiff",
  units = "mm", width = 100, height = 60, dpi = 300
)

# Two-way ANOVA -----------------------------------------------------------

# https://www.sheffield.ac.uk/polopoly_fs/1.536444!/file/MASH_2way_ANOVA_in_R.pdf
# https://dzchilds.github.io/stats-for-bio/two-way-anova-in-r.html

anova2 <- aov(X664.1167 ~ as.factor(Treatment) * as.factor(Time), data = dataTable)

summary(anova2)

res <- TukeyHSD(anova2)

sigTable <- data.frame(res$`as.factor(Treatment):as.factor(Time)`)
sigTable

# CIA:2 weeks-Ctrl:2 weeks 0.0097990 **
# CIA:4 weeks-Ctrl:4 weeks 0.0496138 *
# CIA:6 weeks-Ctrl:6 weeks 0.0173971 *
# CIA:10 weeks-Ctrl:10 weeks 1.0000000

# Save p values
fwrite(sigTable, file = "./result/DMD1659/ind_mz/anova_nad_manuscript.csv", row.names = TRUE)
