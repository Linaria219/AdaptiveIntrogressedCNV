library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(scales)  

file_path <- "path/to/your/cnv_data.txt"

cnv_data <- fread(file_path, header = FALSE, sep = "\t")

# Get basic file information
n_cnvs <- nrow(cnv_data)  # CNV number
n_samples <- ncol(cnv_data) - 4  # Sample number (subtracting the first 4 columns)

# Initialize status matrix
status_matrix <- data.frame(
  row_id = 1:n_cnvs,
  status = rep(-1, n_cnvs),  # Initial status all set to -1
  cnv_type = cnv_data[[4]]  # Fourth column is CNV type
)

# Create result storage data frame
results <- data.frame(
  sample_count = integer(),
  singleton_count = integer(),
  nonsingleton_count = integer(),
  total_detected = integer()
)

# Process each sample
for (sample_idx in 1:n_samples) {
  # Current sample column index (starting from the 5th column)
  current_col <- 4 + sample_idx
  
  # Get CN values for the current sample
  cn_values <- cnv_data[[current_col]]
  
  # Update status matrix
  for (i in 1:n_cnvs) {
    if (cn_values[i] != 2) {  # CN value not equal to 2 indicates CNV detection
      current_status <- status_matrix$status[i]
      
      # Update based on current status
      if (current_status == -1) {
        status_matrix$status[i] <- 0  # First detection
      } else if (current_status == 0) {
        status_matrix$status[i] <- 1  # Second detection
      }
    }
  }
  
  # Statistics for the current status
  singleton_count <- sum(status_matrix$status == 0)
  nonsingleton_count <- sum(status_matrix$status == 1)
  total_detected <- singleton_count + nonsingleton_count
  
  # Record results
  results <- rbind(results, data.frame(
    sample_count = sample_idx,
    singleton_count = singleton_count,
    nonsingleton_count = nonsingleton_count,
    total_detected = total_detected
  ))
  
  # Progress update
  if (sample_idx %% 100 == 0) {
    cat("Processed", sample_idx, "/", n_samples, "samples\n")
  }
}

# Data reshaping to long format for ggplot visualization
results_long <- results %>%
  select(sample_count, singleton_count, nonsingleton_count) %>%
  pivot_longer(cols = c(singleton_count, nonsingleton_count),
               names_to = "cnv_category",
               values_to = "count")

# Set category labels
results_long$cnv_category <- factor(results_long$cnv_category,
                                    levels = c("nonsingleton_count", "singleton_count"),
                                    labels = c("Non-Singleton", "Singleton"))

academic_colors <- c(
  "Non-Singleton" = "#0072B2",  
  "Singleton" = "#D55E00"       
)

# Create stacked area plot
p <- ggplot(results_long, aes(x = sample_count, y = count, fill = cnv_category)) +
  geom_area(alpha = 0.85, linewidth = 0.3, color = "black") +
  scale_fill_manual(values = academic_colors) +
  scale_y_continuous(
    labels = scales::comma_format(),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "CNV Detection Patterns with Increasing Sample Size",
    subtitle = "Accumulation of singleton and non-singleton CNVs across samples",
    x = "Number of Samples",
    y = "Number of CNVs",
    fill = "CNV Category"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", face = "bold"),
    plot.title = element_text(
      hjust = 0.5, 
      size = 20, 
      face = "bold", 
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(
      hjust = 0.5, 
      size = 16, 
      face = "plain",  
      margin = margin(b = 20)
    ),
    axis.title = element_text(
      size = 18, 
      face = "bold"
    ),
    axis.text = element_text(
      size = 16, 
      face = "bold"
    ),
    axis.text.x = element_text(
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(
      margin = margin(r = 5)
    ),
    legend.title = element_text(
      size = 18, 
      face = "bold"
    ),
    legend.text = element_text(
      size = 16, 
      face = "bold"
    ),
    legend.position = "top",
    legend.box = "horizontal",
    legend.spacing = unit(0.5, "cm"),
    legend.key.size = unit(1.5, "cm"),
    panel.grid.major = element_line(
      color = "#E0E0E0",  
      linewidth = 0.5,
      linetype = "solid"
    ),
    panel.grid.minor = element_line(
      color = "#F0F0F0",  
      linewidth = 0.3,
      linetype = "solid"
    ),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(30, 30, 30, 30, "pt"),
    panel.border = element_rect(
      color = "#D0D0D0",  
      fill = NA,
      linewidth = 0.8
    )
  ) +
  guides(fill = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    nrow = 1,
    byrow = TRUE
  ))

print(p)

png_output <- "CNV_Detection_Patterns.png"
png(png_output, width = 4000, height = 3000, res = 300, 
    family = "Helvetica", pointsize = 20)
print(p)
dev.off()

svg_output <- "CNV_Detection_Patterns.svg"
svg(svg_output, width = 20, height = 15, 
    family = "Helvetica", pointsize = 20)
print(p)
dev.off()