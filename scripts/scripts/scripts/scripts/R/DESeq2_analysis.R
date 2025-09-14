# DESeq2_analysis.R - Differential expression analysis
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)
library(RColorBrewer)

# Paths
fc_file <- "results/counts/featurecounts.txt"
meta_file <- "data/sample_info.csv"
out_dir <- "results/deseq2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read featureCounts
fc <- read.delim(fc_file, comment.char="#", stringsAsFactors=FALSE)
counts <- fc[, 7:ncol(fc)]
rownames(counts) <- fc$Geneid
colnames(counts) <- gsub(".sorted.bam$|.bam$","", colnames(counts))

# Read metadata
coldata <- read.csv(meta_file, row.names = 1, stringsAsFactors=FALSE)
coldata <- coldata[colnames(counts), , drop=FALSE]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = coldata, design = ~ condition)
dds <- dds[ rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Results
res <- results(dds, alpha = 0.05)
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), file = file.path(out_dir, "deseq2_results_all.csv"))

# LFC shrink
print(resultsNames(dds))
resLFC <- lfcShrink(dds, contrast = c("condition","treated","control"), type="apeglm")
write.csv(as.data.frame(resLFC), file = file.path(out_dir, "deseq2_results_shrunk.csv"))

# PCA
vsd <- vst(dds, blind = FALSE)
pdf(file.path(out_dir, "PCA_plot.pdf"))
plotPCA(vsd, intgroup = "condition")
dev.off()

# Volcano plot
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)
res_df <- na.omit(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, "yes", "no")
res_df$minusLog10Padj <- -log10(res_df$padj + 1e-300)

p <- ggplot(res_df, aes(x = log2FoldChange, y = minusLog10Padj, color = significant)) +
  geom_point(alpha = 0.6, size=1.2) +
  scale_color_manual(values = c("no" = "grey70", "yes" = "red")) +
  theme_minimal() +
  labs(x = "log2 fold change", y = "-log10 (padj)", title = "Volcano plot")
ggsave(filename = file.path(out_dir, "volcano.png"), plot = p, width = 6, height = 4)

# Heatmap top 50 variable genes
rv <- rowVars(assay(vsd))
topgenes <- head(order(rv, decreasing = TRUE), 50)
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)
pdf(file.path(out_dir, "heatmap_top50.pdf"), width=8, height=10)
pheatmap(mat, annotation_col = coldata[, "condition", drop=FALSE], show_rownames = TRUE)
dev.off()
