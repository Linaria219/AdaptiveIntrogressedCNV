library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(scales)

# Set target region
target_chr <- "chr14"
region_start <- 106031722
region_end <- 106300478
region_length <- region_end - region_start

cat("Target region:", target_chr, ":", 
    format(region_start, big.mark = ","), "-",
    format(region_end, big.mark = ","), "\n")
cat("Region length:", format(region_length, big.mark = ","), "bp\n")

# Extract gene information from GTF file
gtf_file <- "/path/to/your/annotation.gtf"

# Read GTF file
gtf_data <- rtracklayer::import(gtf_file)

# Only keep exon rows (GeneRegionTrack requires exon information)
exon_gtf <- gtf_data[gtf_data$type == "exon"]

region_gr <- GRanges(seqnames = target_chr,
                     ranges = IRanges(start = region_start, end = region_end))

# Filter exons in the target region
exons_in_region <- subsetByOverlaps(exon_gtf, region_gr)

cat("In the target region, found", length(unique(exons_in_region$gene_id)), "genes with", 
    length(exons_in_region), "exons\n")

# Prepare data format required by Gviz
if (length(exons_in_region) > 0) {
  # Convert to data frame
  exon_df <- as.data.frame(exons_in_region)
  
  # Rename columns to meet Gviz requirements
    gene_models <- data.frame(
    chromosome = exon_df$seqnames,
    start = exon_df$start,
    end = exon_df$end,
    width = exon_df$width,
    strand = exon_df$strand,
    feature = "exon",
    gene = ifelse(!is.null(exon_df$gene_name) & !is.na(exon_df$gene_name),
                  exon_df$gene_name, exon_df$gene_id),
    exon = exon_df$exon_id,
    transcript = exon_df$transcript_id,
    symbol = ifelse(!is.null(exon_df$gene_name) & !is.na(exon_df$gene_name),
                    exon_df$gene_name, exon_df$gene_id),
    stringsAsFactors = FALSE
  )
  
  # Remove duplicate rows
  gene_models <- unique(gene_models)
  
  # Print gene statistics
  unique_genes <- unique(gene_models$gene)
  cat("Genes in the region:", paste(unique_genes, collapse = ", "), "\n")
  
} else {
  cat("No exons found in the specified region\n")
  stop("No exons found, cannot proceed with plotting.")
}

# Create Gviz tracks
cat("\nCreating Gviz tracks...\n")

# Set UCSC chromosome name options (Gviz requires chromosomes to start with "chr")
options(ucscChromosomeNames = FALSE)

# Create ideogram track
cat("Creating ideogram track...\n")
ideoTrack <- tryCatch({
  IdeogramTrack(genome = "hg19", chromosome = target_chr)
}, error = function(e) {
  cat("Failed to create IdeogramTrack, using alternative approach...\n")
  # If you are unable to create an IdeogramTrack, use a simple annotated track instead.
  NULL
})

# Create genome axis track
cat("Creating genome axis track...\n")
axisTrack <- GenomeAxisTrack(
  labelPos = "alternating",
  add53 = TRUE,
  add35 = TRUE,
  littleTicks = TRUE,
  fontcolor = "black",
  fontsize = 16,  
  fontface = 1,   
  fontfamily = "sans"  
)

# Create gene region track (core track)
cat("Creating gene region track...\n")
geneTrack <- GeneRegionTrack(
  range = gene_models,
  genome = "hg19",
  chromosome = target_chr,
  name = "Gene structure",
  transcriptAnnotation = "symbol",
  background.title = "darkblue",
  col.title = "white",
  showId = TRUE,
  cex.title = 1.0,  
  cex.group = 0.9,  
  fontcolor.group = "black",
  fontface.group = 1,
  fontfamily.group = "sans",  
  col = NULL,  
  fill = "#1E88E5",  
  col.line = "darkblue",  
  fontcolor = "black",
  fontsize = 16,  
  fontface = 1,
  fontfamily = "sans"  
)

# Create gene annotation track (display gene names)
cat("Creating gene annotation track...\n")
annoTrack <- AnnotationTrack(
  start = sapply(split(gene_models$start, gene_models$gene), min),
  end = sapply(split(gene_models$end, gene_models$gene), max),
  chromosome = target_chr,
  genome = "hg19",
  name = "Gene",
  id = unique(gene_models$gene),
  group = unique(gene_models$gene),
  stacking = "dense",
  background.title = "darkgreen",
  col.title = "white",
  cex.title = 1.0,  
  cex.feature = 0.9,  
  fontcolor.feature = "black",
  fontface.feature = 1,
  fontfamily.feature = "sans",  
  fontcolor = "black",
  fontsize = 16,  
  fontface = 1,
  fontfamily = "sans"  # Use sans-serif font
)

# The list of tracks to be drawn depends on whether there is an IdeogramTrack.
if (!is.null(ideoTrack)) {
  tracks <- list(ideoTrack, axisTrack, geneTrack, annoTrack)
  track_sizes <- c(0.5, 0.5, 2, 0.5)  
} else {
  tracks <- list(axisTrack, geneTrack, annoTrack)
  track_sizes <- c(0.5, 2, 0.5)
}

# Set graphical parameters, use sans-serif font
par(family = "sans")

plotTracks(
  tracks,
  from = region_start,
  to = region_end,
  chromosome = target_chr,
  sizes = track_sizes,
  main = paste("chr14:", format(region_start, big.mark = ","), 
               "-", format(region_end, big.mark = ","), "region"),
  cex.main = 1.5,  
  background.panel = "#FFFEDB",
  background.title = "darkblue",
  col.axis = "black",
  cex.axis = 1.2,  
  fontface.axis = 1,
  fontfamily.axis = "sans"  
)

output_file <- paste0("gviz_gene_region_arial_", target_chr, "_", 
                      format(region_start, scientific = FALSE), "_",
                      format(region_end, scientific = FALSE), ".png")

png(output_file, width = 1800, height = 1000, res = 150)

plotTracks(
  tracks,
  from = region_start,
  to = region_end,
  chromosome = target_chr,
  sizes = track_sizes,
  main = paste("chr14:", format(region_start, big.mark = ","), 
               "-", format(region_end, big.mark = ","), "region"),
  cex.main = 1.5,
  background.panel = "#FFFEDB",
  background.title = "darkblue",
  col.axis = "black",
  cex.axis = 1.2,
  fontface.axis = 1,
  fontfamily.axis = "sans"
)

dev.off()

cat("The PNG image has been saved as:", output_file, "\n")

svg_file <- paste0("gviz_gene_region_arial_", target_chr, "_", 
                   format(region_start, scientific = FALSE), "_",
                   format(region_end, scientific = FALSE), ".svg")

svg(svg_file, width = 16, height = 9)

plotTracks(
  tracks,
  from = region_start,
  to = region_end,
  chromosome = target_chr,
  sizes = track_sizes,
  main = paste("chr14:", format(region_start, big.mark = ","), 
               "-", format(region_end, big.mark = ","), "region"),
  cex.main = 1.5,
  background.panel = "#FFFEDB",
  background.title = "darkblue",
  col.axis = "black",
  cex.axis = 1.2
)

dev.off()

cat("The SVG image has been saved as:", svg_file, "\n")
