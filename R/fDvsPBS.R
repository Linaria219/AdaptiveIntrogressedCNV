library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)
library(cowplot)

file1_path <- "path/to/your/PBSandfD.bed"
file2_path <- "path/to/your/CandidateCNV.bed"

data1 <- read.table(file1_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

if (ncol(data1) != 15) {
  stop("The first file has an incorrect number of columns; 15 columns were expected, but the actual number is incorrect," ncol(data1), "columns".)
}

colnames(data1) <- c("chr", "start", "end", "type", 
                     "V5", "size_value", paste0("V", 7:15))  # 6th column is size_value

data1 <- data1 %>%
  dplyr::select(chr, start, end, type, size_value,
                pbs_pvalue = V10, 
                denisovan_pvalue = V12, 
                neanderthal_pvalue = V14)

data1$id <- paste(data1$chr, data1$start, data1$end, sep = "_")

data2 <- read.table(file2_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

if (ncol(data2) != 3) {
  stop("The second file has an incorrect number of columns; 3 columns were expected, but the actual number is incorrect," ncol(data2), "columns".)
}

colnames(data2) <- c("chr", "start", "end")
data2$id <- paste(data2$chr, data2$start, data2$end, sep = "_")

# Data cleaning function: handle NA, nan, and scientific notation
clean_pvalues <- function(pvec) {
  # Convert to character
  pvec_char <- as.character(pvec)
  
  # Handle special values
  pvec_clean <- rep(NA, length(pvec_char))
  
  for (i in seq_along(pvec_char)) {
    val <- pvec_char[i]
    
    # Check if it is NA or nan
    if (is.na(val) || tolower(val) == "nan") {
      pvec_clean[i] <- NA
    } else {
      # Try to convert to numeric
      num_val <- tryCatch(
        as.numeric(val),
        warning = function(w) NA,
        error = function(e) NA
      )
      
      pvec_clean[i] <- num_val
    }
  }
  
  return(pvec_clean)
}

# Clean up P-value data and size_value
data1$pbs_pvalue_clean <- clean_pvalues(data1$pbs_pvalue)
data1$denisovan_pvalue_clean <- clean_pvalues(data1$denisovan_pvalue)
data1$neanderthal_pvalue_clean <- clean_pvalues(data1$neanderthal_pvalue)
data1$size_value_clean <- clean_pvalues(data1$size_value)

# Check the cleaning results
cat("Valid PBS P Value data:", sum(!is.na(data1$pbs_pvalue_clean)), "/", nrow(data1), "\n")
cat("Valid Denisovan P Value data:", sum(!is.na(data1$denisovan_pvalue_clean)), "/", nrow(data1), "\n")
cat("Valid Neanderthal P Value data:", sum(!is.na(data1$neanderthal_pvalue_clean)), "/", nrow(data1), "\n")
cat("Valid Size value data:", sum(!is.na(data1$size_value_clean)), "/", nrow(data1), "\n")

# Remove invalid P-values (rows where all P-values are NA)
data1_clean <- data1 %>%
  filter(!(is.na(pbs_pvalue_clean) & is.na(denisovan_pvalue_clean) & is.na(neanderthal_pvalue_clean)))

cat("Number of rows after cleaning:", nrow(data1_clean), "\n")

# Define CNV types
cnv_types <- c("DEL", "DUP", "mCNV")

# =====================================================================
# Custom coordinate transformation functions
# =====================================================================

# X-axis transformation function: 1 on the right, 0.005 at 1/3 position, 0.0001 at 2/3 position
x_transform <- function(p) {
  # Make sure p is within [0,1]
  p <- pmax(pmin(p, 1), 0)
  
  # Define key points
  key_points <- c(1, 0.005, 0.0001, 0)
  positions <- c(1, 5/6, 2/3, 0)
  
  # Linear interpolation
  if (p >= 1) return(1)
  if (p <= 0) return(0)
  
  # Find the interval where p lies
  for (i in 1:(length(key_points)-1)) {
    if (p <= key_points[i] && p >= key_points[i+1]) {
      # Linear interpolation
      t <- (p - key_points[i]) / (key_points[i+1] - key_points[i])
      return(positions[i] + t * (positions[i+1] - positions[i]))
    }
  }
  
  # If p is not in any interval, return the closest point
  return(positions[which.min(abs(key_points - p))])
}

# Y-axis transformation function (reversed): 1 at the bottom, 0.05 at 1/3 position, 0.001 at 2/3 position
y_transform <- function(p) {
  # Make sure p is within [0,1]
  p <- pmax(pmin(p, 1), 0)
  
  # Define key points
  key_points <- c(1, 0.05, 0.001, 0)
  positions <- c(0, 1/6, 1/3, 1)
  
  # Linear interpolation
  if (p >= 1) return(0)
  if (p <= 0) return(1)
  
  # Find the interval where p lies
  for (i in 1:(length(key_points)-1)) {
    if (p <= key_points[i] && p >= key_points[i+1]) {
      # Linear interpolation
      t <- (p - key_points[i]) / (key_points[i+1] - key_points[i])
      return(positions[i] + t * (positions[i+1] - positions[i]))
    }
  }
  
  # If p is not in any interval, return the closest point
  return(positions[which.min(abs(key_points - p))])
}

# Create vectorized transformation functions
y_transform_vec <- Vectorize(y_transform)
x_transform_vec <- Vectorize(x_transform)

# =====================================================================
# Create Adaptive introgression candidate classification
# =====================================================================

# Mark Adaptive introgression candidates
cat("Classifying Adaptive introgression candidates...\n")

# First, mark the basic Adaptive candidates (present in file2)
data1_clean$in_adaptive_file <- data1_clean$id %in% data2$id

# Classify each CNV
data1_clean$adaptive_status <- "non-adaptive"

for (i in 1:nrow(data1_clean)) {
  in_file <- data1_clean$in_adaptive_file[i]
  denisovan_p <- data1_clean$denisovan_pvalue_clean[i]
  neanderthal_p <- data1_clean$neanderthal_pvalue_clean[i]
  
  if (in_file) {
    # Check if Denisovan and Neanderthal P-values are both less than 0.05
    denisovan_sig <- !is.na(denisovan_p) && denisovan_p < 0.05
    neanderthal_sig <- !is.na(neanderthal_p) && neanderthal_p < 0.05
    
    if (denisovan_sig && neanderthal_sig) {
      data1_clean$adaptive_status[i] <- "both_archaic"
    } else if (denisovan_sig) {
      data1_clean$adaptive_status[i] <- "denisovan_only"
    } else if (neanderthal_sig) {
      data1_clean$adaptive_status[i] <- "neanderthal_only"
    } else {
      # In file but P-values are not less than 0.05
      data1_clean$adaptive_status[i] <- "in_file_only"
    }
  }
}

# Create size_category: classify based on size_value_clean values
data1_clean$size_category <- cut(
  data1_clean$size_value_clean,
  breaks = c(-Inf, 0.1, 0.5, Inf),
  labels = c("small", "medium", "large"),
  right = FALSE,
  include.lowest = TRUE
)

# Set the point size mapping
size_mapping <- c(
  "small" = 2.0,   # Minimum radius
  "medium" = 2.5,  # Medium radius
  "large" = 3.0    # Maximum radius
)

# Output classification statistics
cat("Adaptive introgression candidate classification:\n")
cat("  non-adaptive:", sum(data1_clean$adaptive_status == "non-adaptive"), "\n")
cat("  in_file_only:", sum(data1_clean$adaptive_status == "in_file_only"), "\n")
cat("  denisovan_only:", sum(data1_clean$adaptive_status == "denisovan_only"), "\n")
cat("  neanderthal_only:", sum(data1_clean$adaptive_status == "neanderthal_only"), "\n")
cat("  both_archaic:", sum(data1_clean$adaptive_status == "both_archaic"), "\n")

cat("\nSize value classification:\n")
cat("  small (<0.1):", sum(data1_clean$size_category == "small", na.rm = TRUE), "\n")
cat("  medium (0.1-0.5):", sum(data1_clean$size_category == "medium", na.rm = TRUE), "\n")
cat("  large (>=0.5):", sum(data1_clean$size_category == "large", na.rm = TRUE), "\n")

# =====================================================================
# Create plotting function
# =====================================================================
create_plot_custom <- function(data, cnv_type, y_var) {
  # Filter specific CNV types
  plot_data <- data %>% filter(type == cnv_type)
  
  # Select the appropriate P-values and labels based on y_var
  if (y_var == "denisovan") {
    plot_data <- plot_data %>% 
      filter(!is.na(denisovan_pvalue_clean) & !is.na(pbs_pvalue_clean))
    y_values <- plot_data$denisovan_pvalue_clean
    y_label <- "Denisovan introgression P value"
    
    # For the Denisovan plot, select the relevant Adaptive candidates
    # Including denisovan_only and both_archaic
    plot_data$point_color <- ifelse(
      plot_data$adaptive_status %in% c("denisovan_only", "both_archaic"),
      "adaptive",
      "non-adaptive"
    )
    
    # Further distinguish both_archaic
    plot_data$is_both_archaic <- plot_data$adaptive_status == "both_archaic"
  } else if (y_var == "neanderthal") {
    plot_data <- plot_data %>% 
      filter(!is.na(neanderthal_pvalue_clean) & !is.na(pbs_pvalue_clean))
    y_values <- plot_data$neanderthal_pvalue_clean
    y_label <- "Neanderthal introgression P value"
    
    # For the Neanderthal plot, select the relevant Adaptive candidates
    # Including neanderthal_only and both_archaic
    plot_data$point_color <- ifelse(
      plot_data$adaptive_status %in% c("neanderthal_only", "both_archaic"),
      "adaptive",
      "non-adaptive"
    )
    
    # Further distinguish both_archaic
    plot_data$is_both_archaic <- plot_data$adaptive_status == "both_archaic"
  }
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = paste("No data for", cnv_type, y_var), 
                      family = "Helvetica", size = 6, face = "bold") +
             theme_void())
  }
  
  # Apply custom transformation
  plot_data$x_transformed <- x_transform_vec(plot_data$pbs_pvalue_clean)
  plot_data$y_transformed <- y_transform_vec(y_values)
  
  # Calculate the transformed positions of the threshold lines
  threshold_y <- y_transform(0.05)
  threshold_x <- x_transform(0.005)
  
  # Assign default values for missing size_category data
  plot_data$size_category <- factor(
    plot_data$size_category,
    levels = c("small", "medium", "large")
  )
  
   p <- ggplot(plot_data, aes(x = x_transformed, y = y_transformed)) +
    geom_point(data = plot_data %>% filter(point_color == "non-adaptive"),
               aes(size = size_category),
               color = "#CCCCCC", alpha = 0.5, shape = 19) +
    
    geom_point(data = plot_data %>% filter(point_color == "adaptive" & !is_both_archaic),
               aes(size = size_category),
               color = "#4E79A7", alpha = 0.5, shape = 19) +
    
    geom_point(data = plot_data %>% filter(is_both_archaic),
               aes(size = size_category),
               color = "#FF5500", alpha = 0.5, shape = 19) +
    
    geom_hline(yintercept = threshold_y, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.5) +
    
    geom_vline(xintercept = threshold_x, linetype = "dashed", color = "blue", alpha = 0.7, linewidth = 0.5) +
    
    scale_x_continuous(
      breaks = c(1, 5/6, 2/3, 0),
      labels = c("1", "0.005", "0.0001", "0"),
      limits = c(0, 1),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      breaks = c(0, 1/6, 1/3, 1),
      labels = c("1", "0.05", "0.001", "0"),
      limits = c(0, 1),
      expand = expansion(mult = 0.02)
    ) +
    
    scale_size_manual(
      values = size_mapping,
      name = "Size value",
      labels = c("small" = "<0.1", "medium" = "0.1-0.5", "large" = "≥0.5")
    ) +
    
    labs(
      title = paste(cnv_type, "-", ifelse(y_var == "denisovan", "Denisovan", "Neanderthal")),
      x = "PBS P value",
      y = y_label
    ) +

    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", family = "Helvetica"),  
      axis.title = element_text(size = 16, face = "bold", family = "Helvetica"),  
      axis.text = element_text(size = 14, face = "bold", family = "Helvetica"),  
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "none"  #
    )
  
  return(p)
}

# =====================================================================
# Create and save plots
# =====================================================================

# Create all 6 plots

# Top row: Denisovan P-value
p1 <- create_plot_custom(data1_clean, "DEL", "denisovan")
p2 <- create_plot_custom(data1_clean, "DUP", "denisovan")
p3 <- create_plot_custom(data1_clean, "mCNV", "denisovan")

# Bottom row: Neanderthal P-value
p4 <- create_plot_custom(data1_clean, "DEL", "neanderthal")
p5 <- create_plot_custom(data1_clean, "DUP", "neanderthal")
p6 <- create_plot_custom(data1_clean, "mCNV", "neanderthal")

# Create a layout matrix
layout_matrix <- rbind(c(1, 2, 3), c(4, 5, 6))

# Use grid.arrange to combine the plots
combined_plot <- grid.arrange(
  p1, p2, p3, p4, p5, p6,
  nrow = 2, ncol = 3,
  layout_matrix = layout_matrix,
  top = textGrob("PBS P value vs Introgression P values across CNV Types", 
                 gp = gpar(fontsize = 22, fontface = "bold", fontfamily = "Helvetica"))  
)

# =====================================================================
create_global_legend <- function() {
  color_legend_data <- data.frame(
    type = c("Not adaptive introgression", "Adaptive introgression", "Adaptive introgressed from both Archaic"),
    color = c("#CCCCCC", "#4E79A7", "#FF5500"),
    x = 1:3,
    y = 1
  )
  
  size_legend_data <- data.frame(
    size_category = factor(c("small", "medium", "large"), 
                           levels = c("small", "medium", "large")),
    label = c("<0.1", "0.1-0.5", "≥0.5"),
    x = 1:3,
    y = 0.5
  )
  
  color_legend <- ggplot(color_legend_data, aes(x = x, y = y, color = type)) +
    geom_point(size = 9, shape = 20) +  
    scale_color_manual(values = c("Not adaptive introgression" = "#CCCCCC", 
                                  "Adaptive introgression" = "#4E79A7",
                                  "Adaptive introgressed from both Archaic" = "#FF5500"),
                       name = "Adaptive Status") +
    theme_void() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12, face = "bold", family = "Helvetica"),  
      legend.title = element_text(size = 14, face = "bold", family = "Helvetica")  
    ) +
    guides(color = guide_legend(override.aes = list(size = 8)))  
  
  size_legend <- ggplot(size_legend_data, aes(x = x, y = y, size = size_category)) +
    geom_point(shape = 19, color = "black") +
    scale_size_manual(values = size_mapping,
                      name = "Vst",
                      labels = c("<0.1", "0.1-0.5", "≥0.5")) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 14, face = "bold", family = "Helvetica"),  
      legend.title = element_text(size = 16, face = "bold", family = "Helvetica")  
    )
  
  legend_combined <- plot_grid(
    get_legend(color_legend),
    get_legend(size_legend),
    nrow = 1,
    rel_widths = c(1, 1)
  )
  
  return(legend_combined)
}

create_threshold_legend <- function() {
  threshold_legend <- data.frame(
    type = c("y=0.05 threshold", "x=0.005 threshold"),
    color = c("red", "blue"),
    linetype = "dashed"
  )
  
  p <- ggplot(threshold_legend, aes(x = 0, xend = 1, y = 0, yend = 0, 
                                    color = type, linetype = linetype)) +
    geom_segment(linewidth = 1.2) +  
    scale_color_manual(values = c("y=0.05 threshold" = "red", 
                                  "x=0.005 threshold" = "blue"),
                       name = "Threshold Lines") +
    scale_linetype_manual(values = c("dashed" = "dashed"), guide = "none") +
    theme_void() +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 14, face = "bold", family = "Helvetica"),  
      legend.title = element_text(size = 16, face = "bold", family = "Helvetica")  
    )
  
  return(get_legend(p))
}

legend_combined <- create_global_legend()
threshold_legend <- create_threshold_legend()

full_legend <- plot_grid(
  legend_combined,
  threshold_legend,
  nrow = 2,
  rel_heights = c(1, 0.5)
)

final_plot <- plot_grid(
  combined_plot,
  full_legend,
  nrow = 2,
  rel_heights = c(4, 0.5)
)

print(final_plot)

png_output <- "PBS_vs_Introgression_CNV_Types.png"
tryCatch({
  png(png_output, width = 3200, height = 2200, res = 150)  
  print(final_plot)
  dev.off()
  cat("PNG file has been saved as:", png_output, "\n")
}, error = function(e) {
  stop("Error: can not save PNG file.")
})

svg_output <- "PBS_vs_Introgression_CNV_Types.svg"
tryCatch({
  svg(svg_output, width = 22, height = 15) 
  print(final_plot)
  dev.off()
  cat("SVG file has been saved as:", svg_output, "\n")
}, error = function(e) {
  stop("Error: can not save SVG file.")
})
