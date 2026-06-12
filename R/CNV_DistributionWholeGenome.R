
library(GenomicRanges)
library(karyoploteR)
library(ggplot2)
library(data.table)

file_path <- "path/to/file.bed"

cnv_data <- fread(file_path, header = FALSE, 
                  col.names = c("chr", "start", "end", "type", "length"))

# Filtering autosomes 1-22
cnv_data <- cnv_data[grepl("^chr[1-9]|^chr1[0-9]|^chr2[0-2]$", chr)]

# Convert to GenomicRanges object
cnv_gr <- makeGRangesFromDataFrame(cnv_data, 
                                   keep.extra.columns = TRUE,
                                   seqnames.field = "chr",
                                   start.field = "start", 
                                   end.field = "end")

create_cnv_density_plot <- function(output_file, width, height, res = 300, file_type = "png") {
  if (file_type == "png") {
    png(output_file, width = width, height = height, res = res, 
        family = "Helvetica")
  } else if (file_type == "svg") {
    svg(output_file, width = width/100, height = height/100, 
        family = "Helvetica")
  } else {
    stop("Unsupported file format")
  }

  par(family = "Helvetica")
  
  # Use karyoploteR for genome-wide visualization
  kp <- plotKaryotype(genome = "hg19", 
                      chromosomes = paste0("chr", 1:22),
                      plot.type = 1, 
                      main = "CNV Distribution across Autosomes")  

  colors <- c("DEL" = "#D55E00", "DUP" = "#0072B2", "mCNV" = "#009E73")
  
  # Add CNV density trajectory
  kpAddBaseNumbers(kp)
  
  kpPlotDensity(kp, data = cnv_gr[cnv_gr$type == "DEL"], 
                col = colors["DEL"], 
                window.size = 3e6,      # 2Mb window
                r0 = 0, r1 = 0.3)       # Bottom 30% area
  
  kpPlotDensity(kp, data = cnv_gr[cnv_gr$type == "DUP"], 
                col = colors["DUP"], 
                window.size = 3e6, 
                r0 = 0.3, r1 = 0.6)   # Middle 30% area
  
  kpPlotDensity(kp, data = cnv_gr[cnv_gr$type == "mCNV"], 
                col = colors["mCNV"], 
                window.size = 3e6, 
                r0 = 0.6, r1 = 0.9)     # Top 30% area
  
  # Add legend
  legend("right", legend = names(colors), fill = colors, 
         title = "CNV Types", bty = "n", inset = c(0.15, 0), xpd = TRUE)
  
  # Close graphics device
  dev.off()
}

# Save PNG format
png_output <- "CNV_density_plot.png"
create_cnv_density_plot(png_output, width = 1200, height = 600, res = 300, file_type = "png")

# Save SVG format
svg_output <- "CNV_density_plot.svg"
create_cnv_density_plot(svg_output, width = 1500, height = 750, file_type = "svg")
