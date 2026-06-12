library(ggplot2)
library(data.table)
library(dplyr)
library(scales)

file_path <- "path/to/your/data

data <- fread(file_path, header = FALSE, sep = "\t", 
              col.names = c("chr", "start", "end", "type", "length"),
              colClasses = list(character = 1:4, numeric = 5))

# Data filtering
filtered_data <- data[type %in% c("DEL", "DUP", "mCNV") & length <= 1000000 & length > 0]

rm(data)
gc()

if (nrow(filtered_data) == 0) {
  stop("No data found for the specified types (DEL, DUP, mCNV) with length <= 1,000,000")
}

# Set x-axis ticks and labels
breaks <- c(10, 100, 1000, 10000, 100000, 1000000)
labels <- c("10 bp", "100 bp", "1 kb", "10 kb", "100 kb", "1 Mb")

colors <- c(
  "DEL" = "#D55E00",    
  "DUP" = "#0072B2",    
  "mCNV" = "#009E73"    
)

alphas <- c(
  "DEL" = 0.7,
  "DUP" = 0.6,
  "mCNV" = 0.4
)

p <- ggplot() +
  geom_density(data = filtered_data[type == "DEL"], 
               aes(x = length, color = "DEL", fill = "DEL"), 
               alpha = alphas["DEL"], size = 0.8, adjust = 0.4) +
  geom_density(data = filtered_data[type == "DUP"], 
               aes(x = length, color = "DUP", fill = "DUP"), 
               alpha = alphas["DUP"], size = 0.8, adjust = 0.4) +
  geom_density(data = filtered_data[type == "mCNV"], 
               aes(x = length, color = "mCNV", fill = "mCNV"), 
               alpha = alphas["mCNV"], size = 0.8, adjust = 0.4) +
  
  scale_x_log10(
    breaks = breaks,
    labels = labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  
  scale_color_manual(
    name = "CNV Type",
    values = colors,
    breaks = names(colors)
  ) +
  scale_fill_manual(
    name = "CNV Type",
    values = colors,
    breaks = names(colors)
  ) +
  
  labs(
    x = "CNV Length",
    y = "Density",
    title = "Length Density Distribution of CNV Types"
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
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      margin = margin(r = 10)
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
    # 面板边框
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

png_output <- "CNV_Types_Length_Density_Distribution.png"
png(png_output, width = 4000, height = 2100, res = 300, 
    family = "Helvetica", pointsize = 20)
print(p)
dev.off()

svg_output <- "CNV_Types_Length_Density_Distribution.svg"
svg(svg_output, width = 20, height = 11, 
    family = "Helvetica", pointsize = 20)
print(p)
dev.off()
