library(ggplot2)
library(data.table)
library(dplyr)
library(scales)

file_path <- "path/to/your/cnv_data.bed"  

process_large_bed_file <- function(file_path) {
  cat("Processing large BED file...\n")
  cnv_data <- fread(file_path, sep = "\t", header = TRUE, 
                    select = 1:5,
                    colClasses = list(character = 1:4, numeric = 5))
  cat("File loaded, total rows:", nrow(cnv_data), "\n")
  cnv_data$VAC <- cnv_data[[5]] - 1
  filtered_data <- cnv_data[cnv_data$VAC > 0, ]
  cat("Valid data rows after filtering:", nrow(filtered_data), "\n")
  return(filtered_data)
}

main_analysis <- function() {
  tryCatch({
    cat("=== Using fast reading mode ===\n")
    cnv_data <- process_large_bed_file(file_path)
    if (nrow(cnv_data) > 1000000) {
      cat("Large dataset detected, applying random sampling...\n")
      set.seed(123)
      sample_size <- min(1000000, nrow(cnv_data))
      cnv_data <- cnv_data[sample(1:nrow(cnv_data), sample_size), ]
    }
  }, error = function(e) {
    stop("File reading failed: ", e$message)
  })
  
  # Statistics for each VAC value
  vac_counts <- table(cnv_data$VAC)
  count_data <- data.frame(
    VAC = as.numeric(names(vac_counts)),
    Count = as.numeric(vac_counts)
  )
  
  # Basic statistical information
  cat("\n=== Statistical Summary ===\n")
  cat("Unique VAC values:", nrow(count_data), "\n")
  cat("Total CNVs:", sum(count_data$Count), "\n")
  
  # Generate unified power scale
  max_vac <- max(count_data$VAC)
  max_count <- max(count_data$Count)
  all_breaks <- 10^(0:max(ceiling(log10(max_vac)), ceiling(log10(max_count))))
  
  # Create smooth curve plot (log-log scale, tick labels as regular numbers)
  p <- ggplot(count_data, aes(x = VAC, y = Count)) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#2E86AB", size = .5, span = 0.1) +
    scale_x_log10(
      breaks = all_breaks,
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = all_breaks,
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(
      x = "Variant allele count",
      y = "Number of CNVs"
    ) +
    theme_minimal() +
    theme(
      text = element_text(family = "Helvetica", face = "bold"),
      plot.title = element_text(
        hjust = 0.5, 
        size = 24, 
        face = "bold", 
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        hjust = 0.5, 
        size = 18, 
        face = "plain",  
        margin = margin(b = 20)
      ),
      axis.title = element_text(
        size = 20, 
        face = "bold"
      ),
      axis.text = element_text(
        size = 16, 
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 0, 
        hjust = 0.5, 
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
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing = unit(0.5, "cm"),
      legend.key.size = unit(1.5, "cm"),
      legend.key = element_rect(fill = NA, color = NA),
      panel.grid.major = element_line(
        color = "#E0E0E0",  
        linewidth = 0.5
      ),
      panel.grid.minor = element_line(
        color = "#F0F0F0",  
        linewidth = 0.3
      ),
      panel.border = element_rect(
        color = "#D0D0D0",  
        fill = NA,
        linewidth = 0.8
      ),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(30, 30, 30, 30, "pt")
    )
  
  print(p)
  png_output <- "CNV_VAC_Type_Curves.png"
  png(png_output, width = 4000, height = 3200, res = 300, 
      family = "Helvetica", pointsize = 20)
  print(p)
  dev.off()
  
  svg_output <- "CNV_VAC_Type_Curves.svg"
  svg(svg_output, width = 20, height = 16, 
      family = "Helvetica", pointsize = 20)
  print(p)
  dev.off()
}

result <- main_analysis()
