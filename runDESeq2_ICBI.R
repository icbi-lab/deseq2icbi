#!/usr/bin/env Rscript
'runDESeq2_ICBI.R

Usage:
  runDESeq2_ICBI.R <sample_sheet> <count_table> --result_dir=<res_dir> --c1=<c1> --c2=<c2> [options]
  runDESeq2_ICBI.R --help

Arguments:
  <sample_sheet>                CSV file with the sample annotations.
  <count_table>                 TSV file with the read counts

Mandatory options:
  --result_dir=<res_dir>        Output directory
  --c1=<c1>                     Contrast level 1 (perturbation). Needs to be contained in condition_col.
  --c2=<c2>                     Contrast level 2 (baseline). Needs to be contained in condition_col.

Optional options:
  --replicate_col=<rep_col>     Column that contains the replicate. If <sample_col> not present,
                                assume the sample name = "{cond_col}_R{rep_col}". This corresponds to
                                the nf-core RNA-seq sample-sheet. [default: replicate]
  --condition_col=<cond_col>    Column in sample annotation that contains the condition [default: group]
  --sample_col=<sample_col>     Column in sample annotation that contains the sample names
                                (needs to match the colnames of the count table). Overrides replicate_col.
  --plot_title=<title>          Title shown above plots. Is built from contrast per default.
  --prefix=<prefix>             Results file prefix. Is built from contrasts per default.
  --fdr_cutoff=<fdr>            False discovery rate for GO analysis and volcano plots [default: 0.1]
  --gtf_file=<gtf>              Path to the GTF file used for featurecounts. If specified, a Biotype QC
                                will be performed.
' -> doc

library("conflicted")
library("docopt")
arguments = docopt(doc, version = "0.1")

print(arguments)

library("DESeq2")
library("IHW")
library("org.Hs.eg.db")
library("ggplot2")
library("pcaExplorer")
library("topGO")
library("writexl")
library("readr")
library("dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
library("EnhancedVolcano")
library("ggpubr")
library("tibble")
conflict_prefer("paste", "base")
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)

#### Set some parameters here

# Input and output
sampleAnnotationCSV <- arguments$sample_sheet
readCountFile <- arguments$count_table
results_dir = arguments$result_dir

# prefix and plot title
prefix <- arguments$prefix
if (is.null(prefix)) {
  prefix = paste0(arguments$c1, "_", arguments$c2)
}
plotTitle <- arguments$plot_title
if (is.null(plotTitle)) {
  plotTitle = paste0(arguments$c1, " vs. ", arguments$c2)
}

# Sample information and contrasts
cond_col = arguments$condition_col
replicate_col = arguments$replicate_col
sample_col = arguments$sample_col
contrast = c(cond_col, arguments$c1, arguments$c2)

# Cutoff
fdr_cutoff = as.numeric(arguments$fdr_cutoff)

# GTF for Biotype QC
gtf_file = arguments$gtf_file

# ### Testdata
# sampleAnnotationCSV = "testdata/sampleTableN.csv"
# readCountFile = "testdata/merged_gene_counts.txt"
# results_dir = "out"
# contrast = c("treatment", "PFK158", "DMSO")
# cond_col = "treatment"
# sample_col = "sample"
# replicate_col = "replicate"
# gtf_file = "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf"

############### Start processing
design_formula <- as.formula(paste0("~", cond_col))
sampleAnno <- read_csv(sampleAnnotationCSV) %>%
  filter(get(cond_col) == arguments$c1 | get(cond_col) == arguments$c2)

# Add sample col based on condition and replicate if sample col is not explicitly specified
if(is.null(sample_col)) {
  sample_col = "sample"
  sampleAnno = sampleAnno %>% mutate(sample=paste0(cond_col, "_R", replicate_col))
}


count_mat <- read_tsv(readCountFile) %>%
  mutate(Geneid= remove_ensg_version(Geneid))
ensg_to_genesymbol = count_mat %>% select(Geneid, gene_name)
count_mat = count_mat %>%
  select(c(Geneid, sampleAnno$sample)) %>%
  column_to_rownames("Geneid")

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = sampleAnno,
                              design = design_formula)


## keep only genes where we have >= 10 reads in total
# keep <- rowSums(counts(dds)) >= 10

## keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep,]

# save filtered count file
write_tsv(counts(dds) %>% as_tibble(rownames = "Geneid"), file.path(results_dir, paste0(prefix, "_detectedGenesRawCounts_min_10_reads_in_one_condition.tsv")))

# run DESeq
dds <- DESeq(dds)

# get normalized counts
nc <- counts(dds, normalized=T)

### IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds, filterFun=ihw, contrast=contrast) %>%
  as_tibble(rownames = "Geneid") %>%
  left_join(ensg_to_genesymbol) %>%
  arrange(pvalue)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)


# Filter for adjusted p-value < 0.1
resIHWsig <- resIHW %>% filter(padj < fdr_cutoff)

###### Perform Biotype QC
if(!is.null(gtf_file)) {
  gtf = rtracklayer::import(gtf_file, feature.type="gene") %>%
    as_tibble() %>%
    mutate(gene_id = remove_ensg_version(gene_id))

  count_before = nrow(resIHW)
  resIHW = resIHW %>% left_join(select(gtf, gene_id, gene_type), by=c("Geneid"="gene_id"))
  stopifnot("Number of genes should be the same after adding biotypes"= count_before == nrow(resIHW))
  resIHWsig = resIHWsig %>% left_join(select(gtf, gene_id, gene_type), by=c("Geneid"="gene_id"))

  biotype_counts = resIHWsig %>% group_by(gene_type) %>% count()

  p = biotype_counts %>%
    ggplot(aes(x=gene_type, y=n)) + geom_bar(stat='identity') + theme_bw() + coord_flip()
  ggsave(file.path(results_dir, paste0(prefix, "_biotype_counts.png")), p)
  write_tsv(biotype_counts, file.path(results_dir, paste0(prefix, "_biotype_counts.tsv")))

}



#### write results to TSV and XLSX files
write_tsv(resIHW, file.path(results_dir, paste0(prefix, "_IHWallGenes.tsv")))
write_xlsx(resIHW, file.path(results_dir, paste0(prefix, "_IHWallGenes.xlsx")))


write_tsv(resIHWsig, file.path(results_dir, paste0(prefix, "_IHWsigGenes.tsv")))
write_xlsx(resIHWsig, file.path(results_dir, paste0(prefix, "_IHWsigGenes.xlsx")))






###### Run TOPGO analysis
# significant genes as DE gene, all genes as background
de_symbols <- resIHWsig$Geneid
bg_symbols <- rownames(dds)[rowSums(counts(dds)) > 0]

lapply(c("BP", "MF"), function(ontology) {
  topgoDE <- topGOtable(de_symbols, bg_symbols,
                        ontology = ontology,
                        mapping = "org.Hs.eg.db",
                        geneID = "ENSEMBL")
  write_tsv(topgoDE, file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".tsv")))
  write_xlsx(topgoDE %>% select(-genes), file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".xlsx")))
})




########### PCA plot
vsd <- vst(dds, blind=FALSE)
dds <- estimateSizeFactors(dds)


p <- plotPCA(vsd, intgroup=c(cond_col)) +
  geom_point() +
  geom_text(vjust = 0,hjust = 0.2, nudge_x = -1, nudge_y = 0.5, aes(label = name)) +
  ggtitle(paste0("PCA: ", plotTitle)) +
  scale_color_brewer(type="qual", palette="Set1") +
  theme_bw()
ggsave(file.path(results_dir, paste0(prefix, "_PCA.png")), plot = p)


########## Volcano plot
png(filename =  file.path(results_dir, paste0(prefix, "_vulcano.png")), width = 650, height = 650, units = "px")
EnhancedVolcano(resIHW,
                lab = resIHW$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 1e-6,
                FCcutoff = 2,
                subtitle = "",
                legendPosition = "right",
                title = plotTitle)
dev.off()

png(filename =  file.path(results_dir, paste0(prefix, "_vulcano_padj.png")), width = 650, height = 650, units = "px")
EnhancedVolcano(resIHW,
                lab = resIHW$gene_name,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = fdr_cutoff,
                FCcutoff = 2,
                subtitle = "",
                legendPosition = "right",
                title = plotTitle)
dev.off()
