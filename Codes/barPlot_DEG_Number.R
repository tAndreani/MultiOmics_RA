# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Barplot on number of DEG genes
# Data: Annotated DESeq2 result

# Load library ------------------------------------------------------------

library(ggplot2)
library(data.table)
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

# List annotated DESeq2 result csv file path
csvPath <- list.files("./result/Annotated_DESeq2_Wald_Test_Result", full.names = TRUE)
csvPath

# Import multiple dataframes and bind into one
input <- map_df(csvPath, read.csv)
View(input)

unique(input$Timepoint)

# Get number of DEG -------------------------------------------------------

# Upregulated DEG
up <- input %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange > log2(1.5)) %>%
  group_by(Timepoint) %>%
  summarise(n = n()) %>%
  mutate(Direction = "Up")
up

# Downregulated DEG
down <- input %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange < log2(1 / 1.5)) %>%
  group_by(Timepoint) %>%
  summarise(n = n()) %>%
  mutate(Direction = "Down")
down

# Merge tables
res.deg <- rbind(up, down)
res.deg

# Tidy timepoints
res.deg$Timepoint <- str_replace_all(res.deg$Timepoint, "_", " ")
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
  scale_y_continuous(breaks = c(0, 50, seq(100, 1000, 100))) +
  labs(
    x = NULL, y = "Number of differentially expressed genes",
    caption = "FDR < 0.05 and |log2FoldChange| > 0.58"
  )

ggsave("./plot/DESeq2/barPlot_DEG_Number_rm_pseudogene.tif",
  plot = p, device = "tiff",
  units = "mm", width = 92, height = 60, dpi = 300
)

