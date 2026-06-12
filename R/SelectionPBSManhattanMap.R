library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(viridis)
library(data.table)

file_path <- "path/to/your/data"  

data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("chr", "pos", "pbs")
data$chr_num <- as.numeric(gsub("chr", "", data$chr))

valid_chr <- 1:22
data <- data[data$chr_num %in% valid_chr, ]
data <- data[order(data$chr_num, data$pos), ]

# Calculate the maximum position for each chromosome
chr_max_pos <- data %>%
  group_by(chr_num) %>%
  summarise(max_pos = max(pos, na.rm = TRUE)) %>%
  arrange(chr_num)

# Calculate the offset between chromosomes
chr_offset <- cumsum(c(0, chr_max_pos$max_pos[-nrow(chr_max_pos)] + 1e7))  # Add 1e7 gaps between chromosomes

# Create position mapping
pos_mapping <- data.frame(
  chr_num = chr_max_pos$chr_num,
  offset = chr_offset
)

# Calculate absolute positions for each SNP
data <- data %>%
  left_join(pos_mapping, by = "chr_num") %>%
  mutate(abs_pos = pos + offset)

# Calculate chromosome center positions (for x-axis labels)
chr_centers <- data %>%
  group_by(chr_num) %>%
  summarise(center = mean(abs_pos))

# Calculate Z-scores for PBS values
data$pbs_zscore <- scale(data$pbs)

# Calculate PBS cutoff corresponding to Z-score = 5.45/5.51, using GWAS common threshold 5e-8 and test threshold 1.89e-8
# Formula: cutoff = mean + 5.45 * sd
pbs_mean <- mean(data$pbs, na.rm = TRUE)
pbs_sd <- sd(data$pbs, na.rm = TRUE)
pbs_cutoff <- pbs_mean + 5.51 * pbs_sd

# Calculate the number of SNPs above the cutoff
above_cutoff <- sum(data$pbs > pbs_cutoff, na.rm = TRUE)
total_snps <- nrow(data)
percent_above_cutoff <- round(above_cutoff / total_snps * 100, 2)

viridis_colors <- viridis(22, option = "D")
selected_colors <- viridis_colors

color_mapping <- setNames(selected_colors, valid_chr)

p <- ggplot(data, aes(x = abs_pos, y = pbs, color = factor(chr_num))) +
  geom_point(alpha = 0.7, size = 1.8) +  
  geom_hline(yintercept = pbs_cutoff, 
             color = "red", 
             linetype = "dashed", 
             linewidth = 0.9,
             alpha = 0.8) + 
  scale_color_manual(values = color_mapping, name = "Chromosome") +
  scale_x_continuous(
    breaks = chr_centers$center,
    labels = paste0("chr", chr_centers$chr_num)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))  
  ) +
  labs(
    title = "Genome-wide PBS Manhattan Plot",
    subtitle = paste0("PBS cutoff: ", round(pbs_cutoff, 4), 
                      " (", above_cutoff, " of ", total_snps, 
                      " SNPs, ", percent_above_cutoff, "%)"),
    x = "Chromosome",
    y = "PBS value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Helvetica"),
    plot.subtitle = element_text(hjust = 0.5, size = 16, face = "plain", family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold", family = "Helvetica"),
    axis.text.y = element_text(size = 14, face = "bold", family = "Helvetica"),
    axis.title = element_text(size = 16, face = "bold", family = "Helvetica"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none",  
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 12, face = "bold", family = "Helvetica", color = "gray50"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

print(p)

png_output <- "PBS_Manhattan_Plot_with_Cutoff.png"
tryCatch({
  png(png_output, width = 3200, height = 2000, res = 200)
  print(p)
  dev.off()
  cat("PNG saved successfully\n")
}, error = function(e) {
  cat("Unable to save PNG file, error message:", e$message, "\n")
})

svg_output <- "PBS_Manhattan_Plot_with_Cutoff.svg"
tryCatch({
  svg(svg_output, width = 24, height = 15)
  print(p)
  dev.off()
  cat("SVG saved successfully\n")
}, error = function(e) {
  cat("Unable to save SVG file, error message:", e$message, "\n")
})