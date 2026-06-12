library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)  

file_path <- "path/to/your/area_data.txt"  

# Set parameters for the specified region and sample selection
region_start <- 105981722  # Starting coordinates
region_end <- 106350478    # Ending coordinates
sample_count <- 30         # Randomly select 30 samples
set.seed(123)              # Set random seed for reproducibility

# Specified additional coordinate points
custom_coord1 <- 106031722  # First specified coordinate
custom_coord2 <- 106300478  # Second specified coordinate

# Check if specified coordinates are within the region
if (custom_coord1 < region_start || custom_coord1 > region_end) {
  warning(paste("Custom coordinate 1 (", custom_coord1, ") is not within the specified region (", region_start, "-", region_end, "), using default midpoint"))
  custom_coord1 <- (region_start + region_end) / 2
}

if (custom_coord2 < region_start || custom_coord2 > region_end) {
  warning(paste("Custom coordinate 2 (", custom_coord2, ") is not within the specified region (", region_start, "-", region_end, "), using default midpoint"))
  custom_coord2 <- (region_start + region_end) / 2
}

data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Randomly select sample columns (starting from the 4th column)
sample_col_indices <- 4:ncol(data)
selected_indices <- sample(sample_col_indices, min(sample_count, length(sample_col_indices)))
selected_sample_names <- colnames(data)[selected_indices]

cat("Selected sample count:", length(selected_sample_names), "\n")
cat("First 5 selected samples:", head(selected_sample_names, 5), "\n")

# Extract the required columns: coordinates and selected samples
selected_data <- data[, c(1:3, selected_indices)]
colnames(selected_data)[1:3] <- c("chr", "start", "end")

# Filter data within the specified region
cat("Filtering region:", region_start, "-", region_end, "\n")
region_data <- selected_data %>% 
  filter(start >= region_start & end <= region_end)

# Check if there is any data
if (nrow(region_data) == 0) {
  stop("No data found within the specified region!")
}

cat("Data rows within the region:", nrow(region_data), "\n")

# Debugging: Check coordinate step size
if (nrow(region_data) > 1) {
  steps <- diff(region_data$start)
  cat("Coordinate step size statistics:\n")
  cat("  Step size range:", min(steps), "-", max(steps), "\n")
  cat("  Average step size:", mean(steps), "\n")
  cat("  Step size equal to 500:", all(steps == 500), "\n")
}

# Calculate the expected number of windows
num_windows_expected <- ceiling((region_end - region_start + 1) / 500)
cat("\nExpected window count:", num_windows_expected, "\n")
cat("Actual window count:", nrow(region_data), "\n")

# Check for missing windows
if (nrow(region_data) < num_windows_expected) {
  cat("Warning: Original data has missing windows!\n")
  cat("Missing window count:", num_windows_expected - nrow(region_data), "\n")
}

plot_data_list <- list() # Create a list to store data for each sample

for (i in 1:length(selected_sample_names)) {
  sample_name <- selected_sample_names[i]
  
  sample_data <- region_data[, c("chr", "start", "end", sample_name)]
  colnames(sample_data)[4] <- "CN_value"
  
  sample_data$CN_value <- as.integer(sample_data$CN_value)
  
  sample_data$CN_category <- ifelse(
    sample_data$CN_value >= 10,
    "10+",
    as.character(sample_data$CN_value)
  )
  
  sample_data$Sample <- sample_name
  sample_data$Sample_Index <- i
  
  plot_data_list[[i]] <- sample_data
}

# Merge all sample data
plot_data <- do.call(rbind, plot_data_list)

cn_colors <- c(
  "0" = "#F8F8F8",  
  "1" = "#A0A0A0",  
  "2" = "#000000",  
  "3" = "#00008B",  
  "4" = "#4169E1",  
  "5" = "#1E90FF",  
  "6" = "#20B2AA",  
  "7" = "#32CD32",  
  "8" = "#FFD700",  
  "9" = "#FF69B4",  
  "10+" = "#FF0000" 
)

# Ensure that CN_category is an ordered factor.
cn_levels <- c(as.character(0:9), "10+")
plot_data$CN_category <- factor(plot_data$CN_category, levels = cn_levels)

# Check the actual CN values that appear in the data
cat("\nCN value categories present in the data:\n")
print(table(plot_data$CN_category))

# Use only the colors for CN values that actually appear in the data
present_cn_levels <- levels(plot_data$CN_category)[levels(plot_data$CN_category) %in% unique(plot_data$CN_category)]
present_cn_colors <- cn_colors[names(cn_colors) %in% present_cn_levels]

# Prepare data for rectangle plotting
plot_data$y_min <- plot_data$Sample_Index - 0.4
plot_data$y_max <- plot_data$Sample_Index + 0.4

# Determine x-axis tick positions: start, custom coordinate 1, midpoint, custom coordinate 2, end
x_breaks <- c(region_start, 
              custom_coord1,
              (region_start + region_end) / 2,
              custom_coord2,
              region_end)
x_labels <- format(x_breaks, big.mark = ",", scientific = FALSE)

# Define the function to create the plot
create_cn_plot <- function(font_family = "Arial") {
  p <- ggplot(plot_data) +
    geom_rect(aes(xmin = start, xmax = end, 
                  ymin = y_min, ymax = y_max, 
                  fill = CN_category)) +
    
    scale_fill_manual(values = present_cn_colors, name = "Copy Number", drop = FALSE) +
    
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = 0.01)
    ) +
    
    scale_y_continuous(
      breaks = unique(plot_data$Sample_Index),
      labels = unique(plot_data$Sample),
      expand = expansion(mult = 0.05)
    ) +
    
    labs(
      title = paste("CN Value Distribution Heatmap (", sample_count, "Random Samples)"),
      subtitle = paste("Region: chr14:", format(region_start, big.mark = ","), "-", 
                       format(region_end, big.mark = ","), " | Window Size: 500bp"),
      x = "Chromosome Position (hg19)",
      y = "Samples"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, 
                                margin = margin(b = 10), family = font_family),
      plot.subtitle = element_text(size = 18, face = "bold", hjust = 0.5, 
                                   margin = margin(b = 15), family = font_family),
      axis.title = element_text(size = 18, face = "bold", family = font_family),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.5, 
                                 vjust = 0.5, family = font_family),
      axis.text.y = element_text(size = 14, face = "bold", angle = 0, hjust = 1, family = font_family),
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold", family = font_family),
      legend.text = element_text(size = 14, face = "bold", family = font_family),
      legend.key.size = unit(1.4, "cm"),
      legend.key.height = unit(1.0, "cm"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
  
  return(p)
}

png_output <- paste0("CN_heatmap_", sample_count, "samples_chr14_", 
                     format(region_start, scientific = FALSE), "_",
                     format(region_end, scientific = FALSE), ".png")

tryCatch({
  png(png_output, width = 2800, height = 2000, res = 150)  
  p_arial_png <- create_cn_plot(font_family = "Arial")
  print(p_arial_png)
  dev.off()
  cat("The PNG image has been saved as:", png_output, "\n")
}, error = function(e) {
  stop("Error: Unable to save PNG file with Arial font.")
})

svg_output <- paste0("CN_heatmap_", sample_count, "samples_chr14_", 
                     format(region_start, scientific = FALSE), "_",
                     format(region_end, scientific = FALSE), ".svg")

tryCatch({
  svg(svg_output, width = 24, height = 17)  
  p_arial_svg <- create_cn_plot(font_family = "Arial")
  print(p_arial_svg)
  dev.off()
  cat("The SVG image has been saved as:", svg_output, "\n")
}, error = function(e) {
  stop("Error: Unable to save SVG file with Arial font.")
})