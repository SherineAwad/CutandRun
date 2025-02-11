library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(ReactomePA)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)


#txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


args <- commandArgs(trailingOnly = TRUE)
peakfile <- args[1]
sample <- args[2] 
peaks=GenomicRanges::GRangesList(ZF=readPeakFile(peakfile)) #To pull the merged version and ignore the rest 
seqlevels(peaks) <- paste0("chr", seqlevels(peaks))


peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)


genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
comp  <- compareCluster(
  geneCluster   = genes,
  fun           = "enrichKEGG",
  organism      = "mmu",
  pvalueCutoff  = 1,
  pAdjustMethod = "BH")
   

comp_df <- as.data.frame(comp)
outfile = paste(sample, "compareClusterResult.csv", sep="_")



write.csv(comp_df, outfile, row.names = FALSE)

if(FALSE)
{
#Remove the  - Mus musculus (house mouse suffix 
comp@compareClusterResult$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", comp@compareClusterResult$Description)
#Plot selected pathways
enriched_results <- comp@compareClusterResult
selected_pathway_names <- c("Hippo signaling pathway", "Cell cycle", "MAPK signaling pathway", "Wnt signaling pathway", "AMPK signaling pathway", "TGF-beta signaling pathway", "Phospholipase D signaling pathway", "Polycomb repressive complex", "Notch signaling pathway","PI3K-Akt signaling pathway", "Glioma", "Signaling pathways regulating pluripotency of stem cells") 
selected_pathways <- enriched_results[enriched_results$Description %in% selected_pathway_names, ]

figure_name = paste("ZF", "KEGGSelectedpathways.pdf", sep="_")
pdf(file =figure_name, height=10)
#dotplot(comp, showCategory = selected_pathway_names)
dotplot(comp,  x = "Count", color = "p.adjust",showCategory = selected_pathway_names)
dev.off()


figure_name = paste("ZF", "KEGGpathways.pdf", sep="_")
pdf(file =figure_name)
dotplot(comp , showCategory = 10, title = "KEGG Pathway Enrichment Analysis") + theme(axis.text.x = element_text(size = 8)) 
dev.off()
}


