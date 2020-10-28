library("DESeq2")
library("IHW")
library("org.Hs.eg.db")

library("ggplot2")
library("pcaExplorer")
library("topGO")
library("xlsx")

library("dplyr")
library("EnhancedVolcano")
library(ggpubr)


#### Set some parameters here

# Sample annotation file in TSV format
sampleAnnotationTSV <- "./sampleTableN.tsv"

# read counts file
readCountFile <- "./merged_gene_counts.txt"
  
# Result files prefix
prefix <- "WT_vs_Overexpressed"

# Plot title prefix
plotTitle <- "WT versus Overexpressed"

# set contrast: overexpression or treatment
contrast <- "overexpression"

# set design formula
design_formula <- as.formula(paste0("~",contrast))

#### Done with settings


############### Start processing
sampleAnno <- read.table(sampleAnnotationTSV, row.names = "sample", comment.char = '#', header = TRUE)
coldata <- sampleAnno


cts <- data.matrix(read.table(readCountFile,sep="\t",row.names="Geneid",comment.char = '#', header = TRUE))
counts <- cts[,-c(1)]
counts <- counts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = design_formula)


## keep only genes where we have >= 10 reads in total
# keep <- rowSums(counts(dds)) >= 10

## keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds$treatment))) >= 10
dds <- dds[keep,]

# save filtered count file
write.table(counts(dds),file = paste0(prefix, "_detectedGenesRawCounts_min_10_reads_in_one_condition.tsv"), sep = "\t", quote = FALSE)

# run DESeq
dds <- DESeq(dds)

# get normalized counts
nc <- counts(dds, normalized=T)

### IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)

# order by p-value
resIHWordered <- resIHW[order(resIHW$pvalue),]

# add HGNC Symbols as gene names
ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(resIHWordered),columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
# remove duplicate annotation (main Symbol will be kept)
ann <- ann[!duplicated(ann$ENSEMBL),]

# cbind to DESeq result
resIHWordered.annotated <- cbind(resIHWordered, ann$SYMBOL, ann$GENENAME)

# make dataframe to export as xls
resIHWordered.annotated.df <- cbind(rownames(resIHWordered.annotated), as.data.frame(resIHWordered.annotated))
colnames(resIHWordered.annotated.df) <- c("ENSGENE_ID",colnames(resIHWordered.annotated.df[2:(length(resIHWordered.annotated.df)-2)]),"SYMBOL","DESCRIPTION")
# use ENSGENE_ID if no HGNC symbol available
resIHWordered.annotated.df$SYMBOL <- ifelse(is.na(resIHWordered.annotated.df$SYMBOL), paste(resIHWordered.annotated.df$ENSGENE_ID), resIHWordered.annotated.df$SYMBOL)

# write TSV and XLSX files
write.table(resIHWordered.annotated.df,file = paste0(prefix, "_IHWallGenes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(resIHWordered.annotated.df,file = paste0(prefix, "_IHWallGenes.xlsx"), row.names = FALSE)

# Filter for adjusted p-value < 0.1
resIHWsig <- subset(resIHWordered, padj < 0.1)

# annotate the filtered data as well (TODO make function)
ann <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(resIHWsig),columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
ann <- ann[!duplicated(ann$ENSEMBL),]
resIHWsig.annotated <- cbind(resIHWsig, ann$SYMBOL, ann$GENENAME)
resIHWsig.annotated.df <- cbind(rownames(resIHWsig.annotated), as.data.frame(resIHWsig.annotated))
colnames(resIHWsig.annotated.df) <- c("ENSGENE_ID",colnames(resIHWsig.annotated.df[2:(length(resIHWsig.annotated.df)-2)]),"SYMBOL","DESCRIPTION")
resIHWsig.annotated.df$SYMBOL <- ifelse(is.na(resIHWsig.annotated.df$SYMBOL), paste(resIHWsig.annotated.df$ENSGENE_ID), resIHWsig.annotated.df$SYMBOL)

# write filtered TSV and XLSX files
write.table(resIHWsig.annotated.df,file = paste0(prefix, "_IHWsigGenes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.xlsx(resIHWsig.annotated.df, file = paste0(prefix, "_IHWsigGenes.xlsx"), row.names = FALSE)

###### Run TOPGO analysis
# get gene names
genes <- rownames(resIHWsig)

# significant genes as DE gene, all genes as background
de_symbols <- as.character(genes)
bg_symbols <- rownames(dds)[rowSums(counts(dds)) > 0]

topgoDE <- topGOtable(de_symbols, bg_symbols,
                      ontology = "BP",
                      mapping = "org.Hs.eg.db",
                      geneID = "ENSEMBL")
write.table(topgoDE, file = paste0(prefix, "_topGO_IHWsig.tsv"), sep = "\t", quote = F)
write.xlsx(topgoDE, file = paste0(prefix, "_topGO_IHWsig.xlsx"))



########### PCA plot 
vsd <- vst(dds, blind=FALSE)
dds <- estimateSizeFactors(dds)


p <- plotPCA(vsd, intgroup=c(contrast))
p <- p + geom_point() + geom_text(vjust = 0,hjust = 0.2, nudge_x = -1, nudge_y = 0.5, aes(label = name)) + ggtitle(paste0("PCA: ", plotTitle))
ggsave(paste0(prefix, "_PCA.png"), plot = p)


########## Volcano plot
png(filename = paste0(prefix, "_vulcano.png"), width = 650, height = 650, units = "px")
EnhancedVolcano(resIHWordered.annotated.df,
                lab = resIHWordered.annotated.df$SYMBOL,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 10e-5,
                FCcutoff = 2,
                subtitle = "",
                legendPosition = "right",
                title = plotTitle)
dev.off()

png(filename = paste0(prefix, "_vulcano_padj.png"), width = 650, height = 650, units = "px")
EnhancedVolcano(resIHWordered.annotated.df,
                lab = resIHWordered.annotated.df$SYMBOL,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 2,
                subtitle = "",
                legendPosition = "right",
                title = plotTitle)
dev.off()
