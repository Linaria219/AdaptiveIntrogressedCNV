library(ggplot2)
library(dplyr)
library(viridis)
library(data.table)
library(scales)

file_path <- "path/to/your/data"  

data <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("chr", "start", "end", "fd", "col5", "status")

data$chr_num <- as.numeric(gsub("chr", "", data$chr))

valid_chr <- 1:22
data <- data[data$chr_num %in% valid_chr, ]

# Calculate the center coordinates of the region (the average of the starting and ending coordinates).
data$center <- (data$start + data$end) / 2

# Sort by chromosome and center coordinates
data <- data[order(data$chr_num, data$center), ]

# Calculate the maximum center coordinate for each chromosome
chr_max_pos <- data %>%
  group_by(chr_num) %>%
  summarise(max_pos = max(center, na.rm = TRUE)) %>%
  arrange(chr_num)

# Calculate the offset between chromosomes
chr_offset <- cumsum(c(0, chr_max_pos$max_pos[-nrow(chr_max_pos)] + 5e6))  # Add 5Mb gaps between chromosomes

# Create position mapping
pos_mapping <- data.frame(
  chr_num = chr_max_pos$chr_num,
  offset = chr_offset
)

# Calculate absolute positions for each region
data <- data %>%
  left_join(pos_mapping, by = "chr_num") %>%
  mutate(abs_pos = center + offset)

# Calculate chromosome center positions (for x-axis labels)
chr_centers <- data %>%
  group_by(chr_num) %>%
  summarise(center = mean(abs_pos))

# Calculate Z-scores for fD values
data$fd_zscore <- scale(data$fd)

# Calculate fD cutoff corresponding to Z-score = 4.96
# Formula: cutoff = mean + 4.96 * sd
fd_mean <- mean(data$fd, na.rm = TRUE)
fd_sd <- sd(data$fd, na.rm = TRUE)
#fd_cutoff <- fd_mean + 4.96 * fd_sd
fd_cutoff <- 0.259
  
# Calculate the number of regions above the cutoff
above_cutoff <- sum(data$fd > fd_cutoff, na.rm = TRUE)
total_regions <- nrow(data)
percent_above_cutoff <- round(above_cutoff / total_regions * 100, 2)

viridis_colors_optionC <- viridis(22, option = "C")  # plasma
selected_colors <- viridis_colors_optionC

color_mapping <- setNames(selected_colors, valid_chr)

p <- ggplot(data, aes(x = abs_pos, y = fd, color = factor(chr_num))) +
  geom_point(alpha = 0.7, size = 1.8) + 
  geom_hline(yintercept = fd_cutoff, 
             color = "red", 
             linetype = "dashed", 
             linewidth = 0.9,
             alpha = 0.8) +
  scale_color_manual(values = color_mapping, name = "Chromosome") +
  scale_x_continuous(
    breaks = chr_centers$center,
    labels = paste0("chr", chr_centers$chr_num)
  ) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(
    title = "Genome-wide fD Manhattan Plot",
    x = "Chromosome",
    y = "fD value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", family = "Helvetica"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold", family = "Helvetica"),
    axis.text.y = element_text(size = 16, face = "bold", family = "Helvetica"),
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

png_output <- "fD_Manhattan_Plot_with_Cutoff.nean.png"
tryCatch({
  png(png_output, width = 3000, height = 1800, res = 150)
  print(p)
  dev.off()
  cat("PNG saved successfully\n")
}, error = function(e) {
  cat("Unable to save PNG file, error message:", e$message, "\n")
})

svg_output <- "fD_Manhattan_Plot_with_Cutoff.nean.svg"
tryCatch({
  svg(svg_output, width = 20, height = 12)
  print(p)
  dev.off()
  cat("SVG saved successfully\n")
}, error = function(e) {
  cat("Unable to save SVG file, error message:", e$message, "\n")
})