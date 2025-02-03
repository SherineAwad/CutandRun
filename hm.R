
require(EnrichedHeatmap)
require(rtracklayer)
require(circlize)
require(data.table)

args <- commandArgs(trailingOnly = TRUE)
peak <- args[1]
bw <- args[2] 
split_string <- strsplit(bw, ".bw")[[1]]
mysample <- split_string[1]

targets <- makeGRangesFromDataFrame(df = fread(peak, header = FALSE, data.table = FALSE),seqnames.field = "V1", start.field = "V2", end.field = "V3")
targets <- head(targets, 10000)
ExtendSize <- 1000
targets.extended  <- resize(targets, fix = "center", width = ExtendSize*2)
BigWig <- rtracklayer::import(bw, 
                              format = "BigWig", selection = BigWigSelection(targets.extended))
normMatrix <- normalizeToMatrix(signal = BigWig, target = resize(targets, fix = "center", width = 1),  background = 0, keep = c(0, 0.99),  target_ratio = 0, mean_mode = "w0",  value_column = "score", extend = ExtendSize)

col_fun <- colorRamp2(c(0, 30), c("darkblue", "darkgoldenrod1")) 


title = paste(mysample, " Heatmap", sep =" ") 
EH <- EnrichedHeatmap( mat = normMatrix,pos_line = FALSE, col = col_fun,  column_title = title,
 column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
 use_raster = TRUE, raster_quality = 10, raster_device = "png",
rect_gp = gpar(col = "transparent"), heatmap_legend_param = list(
                         legend_direction = "horizontal",
                         title = "normalized counts"),

 border = FALSE,  top_annotation = HeatmapAnnotation(
  enriched = anno_enriched(
                           gp = gpar(col = "black", lty = 1, lwd=2),
                           col="black")))


figure_name = paste(mysample, "QC.pdf", sep="_")

pdf(figure_name)
draw(EH,   heatmap_legend_side = "bottom",  annotation_legend_side = "bottom",  padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()


