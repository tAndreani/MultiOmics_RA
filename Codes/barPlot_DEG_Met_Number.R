# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Barplot on number of DEG metabolites
# Data: FC and q-value results

# Load library ------------------------------------------------------------

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(plyr)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  axis.title.x = element_blank(),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Load data ---------------------------------------------------------------

input <- fread("./data/IPA_Spotfire_Input/Metabolon_log2FC_qValue_notRounded_Spotfire_20210708.csv")
View(input)

unique(input$Timepoint)

# Get number of DEG -------------------------------------------------------

# Upregulated DEG
up <- input %>%
  filter(FDR < 0.05 & log2FoldChange > log2(1.5)) %>%
  group_by(Timepoint) %>%
  summarise(n = n()) %>%
  mutate(Direction = "Up")
up

# Downregulated DEG
down <- input %>%
  filter(FDR < 0.05 & log2FoldChange < log2(1 / 1.5)) %>%
  group_by(Timepoint) %>%
  summarise(n = n()) %>%
  mutate(Direction = "Down")
down

# Merge tables
res.deg <- rbind(up, down)
res.deg

# Tidy timepoints
res.deg$Timepoint <- factor(res.deg$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Order directions
res.deg$Direction <- factor(res.deg$Direction, levels = c("Up", "Down"))

# Barplot -----------------------------------------------------------------

# Calculate the cumulative sum of n for each timepoint
res.cumsum <- ddply(res.deg, "Timepoint",
                    transform,
                    label_ypos = cumsum(n)
)
head(res.cumsum)

# Stacked barplot
p <- res.cumsum %>%
  ggplot(aes(x = Timepoint, y = n, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(breaks= seq(0, 250, 50)) +
  labs(x = NULL, y = "Number of significant metabolites", 
       caption = "FDR < 0.05 and |log2FoldChange| > 0.58")

ggsave("./plot/barPlot/barPlot_Sig_Metabolite_Number.tif", plot = p, device = "tiff",
       units = "mm", width = 92, height = 60, dpi = 300)
