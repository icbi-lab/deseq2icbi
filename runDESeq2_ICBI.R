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
  --nfcore                      Indicate that the input samplesheet is from the nf-core RNA-seq ppipeline.
                                Will merge entries from the same sample and infer the sample_id from `group` and `replicate`.
                                If set, this option overrides `sample_col`.
  --condition_col=<cond_col>    Column in sample annotation that contains the condition [default: group]
  --sample_col=<sample_col>     Column in sample annotation that contains the sample names
                                (needs to match the colnames of the count table). [default: sample]
  --paired_grp=<paired_grp>     Column that conatins the name of the paired samples, when dealing with
                                paired data.
  --covariate_formula=<formula> Formula to model additional covariates (need to be columns in the samplesheet)
                                that will be appended to the formula built from `condition_col`.
                                E.g. `+ age + sex`. Per default, no covariates are modelled.
  --plot_title=<title>          Title shown above plots. Is built from contrast per default.
  --prefix=<prefix>             Results file prefix. Is built from contrasts per default.
  --fdr_cutoff=<fdr>            False discovery rate for GO analysis and volcano plots [default: 0.1]
  --fc_cutoff=<log2 fc cutoff>  Fold change (log2) cutoff for volcano plots [default: 1]
  --gtf_file=<gtf>              Path to the GTF file used for featurecounts. If specified, a Biotype QC
                                will be performed.
  --gene_id_type=<id_type>      Type of the identifier in the `gene_id` column compatible with AnnotationDbi [default: ENSEMBL]
  --n_cpus=<n_cpus>             Number of cores to use for DESeq2 [default: 1]
  --skip_gsea                   Skip Gene-Set-Enrichment-Analysis step
' -> doc

library("conflicted")
library("docopt")
arguments = docopt(doc, version = "0.1")

print(arguments)

library("BiocParallel")
library("DESeq2")
library("IHW")
library("org.Hs.eg.db")
library("ggplot2")
library("pcaExplorer")
library("topGO")
library("clusterProfiler")
library("ReactomePA")
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
conflict_prefer("rename", "dplyr")
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)

#### Get parameters from docopt

# Input and output
sampleAnnotationCSV <- arguments$sample_sheet
readCountFile <- arguments$count_table
results_dir = arguments$result_dir
paired_grp <- arguments$paired_grp

# prefix and plot title
prefix <- arguments$prefix
plot_title <- arguments$plot_title

# Sample information and contrasts
nfcore = arguments$nfcore
cond_col = arguments$condition_col
sample_col = arguments$sample_col
contrast = c(cond_col, arguments$c1, arguments$c2)
gene_id_type = arguments$gene_id_type
covariate_formula = arguments$covariate_formula

# Cutoff
fdr_cutoff = as.numeric(arguments$fdr_cutoff)
fc_cutoff = as.numeric(arguments$fc_cutoff)

# GTF for Biotype QC
gtf_file = arguments$gtf_file

# Other
n_cpus = as.numeric(arguments$n_cpus)
skip_gsea = arguments$skip_gsea

# # Testdata
# ## Example1
# sampleAnnotationCSV = "testdata/example1/sampleTableN.csv"
# readCountFile = "testdata/example1/merged_gene_counts.txt"
# results_dir = "out"
# paired_grp = NULL
# prefix = "example1"
# plot_title = NULL
# nfcore=FALSE
# cond_col = "treatment"
# sample_col = "sample"
# contrast = c("treatment", "PFK158", "DMSO")
# gene_id_type = "ENSEMBL"
# covariate_formula = ""
# fdr_cutoff = 0.1
# fc_cutoff = 1
# gtf_file = "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf"
# n_cpus = 1
# skip_gsea = FALSE

# ## example_nfcore
# sampleAnnotationCSV = "testdata/example_nfcore/rnaseq_samplesheet.csv"
# readCountFile = "testdata/example_nfcore/salmon.merged.gene_counts.subset.tsv"
# results_dir = "/home/sturm/Downloads/tmp_out"
# nfcore = TRUE
# paired_grp = "donor"
# prefix = NULL
# plot_title = NULL
# cond_col = "group"
# sample_col = NULL
# contrast = c("group", "grpA", "grpB")
# fdr_cutoff = 0.1
# fc_cutoff = 1
# gtf_file = "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf"

############### Sanitize parameters and read input data
register(MulticoreParam(workers = n_cpus))

if (is.null(plot_title)) {
  plot_title = paste0(contrast[[2]], " vs. ", contrast[[3]])
}
if (is.null(prefix)) {
  prefix = paste0(contrast[[2]], "_", contrast[[3]])
}
if (is.null(covariate_formula)) {
  covariate_formula = ""
}

if(is.null(paired_grp)) {
  design_formula <- as.formula(paste0("~", cond_col, covariate_formula))
} else {
  design_formula <- as.formula(paste0("~", paired_grp , " + ", cond_col, covariate_formula))
}

sampleAnno <- read_csv(sampleAnnotationCSV) %>%
  filter(get(cond_col) %in% contrast[2:3])

# Add sample col based on condition and replicate if sample col is not explicitly specified
# and make samplesheet distinct (in case the 'merge replicates' functionality was used).
if(nfcore) {
  sample_col = "sample"
  sampleAnno = sampleAnno %>%
    select(-fastq_1, -fastq_2) %>%
    mutate(sample=paste0(sampleAnno[[cond_col]], "_R", sampleAnno[[replicate_col]])) %>%
    distinct()
}

count_mat <- read_tsv(readCountFile)
if (gene_id_type == "ENSEMBL") {
  count_mat = count_mat %>% mutate(gene_id= remove_ensg_version(gene_id))
}

ensg_to_genesymbol = count_mat %>% select(gene_id, gene_name)
ensg_to_desc = AnnotationDbi::select(org.Hs.eg.db, count_mat$gene_id %>% unique(), keytype = gene_id_type, columns = c("GENENAME")) %>%
  distinct(across(!!gene_id_type), .keep_all = TRUE)

count_mat = count_mat %>%
  select(c(gene_id, sampleAnno[[sample_col]])) %>%
  column_to_rownames("gene_id") %>%
  round() # salmon does not necessarily contain integers


################# Start processing
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = sampleAnno,
                              design = design_formula)


## keep only genes where we have >= 10 reads in total
# keep <- rowSums(counts(dds)) >= 10

## keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep,]

# save filtered count file
write_tsv(counts(dds) %>% as_tibble(rownames = "gene_id"), file.path(results_dir, paste0(prefix, "_detectedGenesRawCounts_min_10_reads_in_one_condition.tsv")))

# save normalized filtered count file
dds <- estimateSizeFactors(dds)
write_tsv(counts(dds, normalized=TRUE) %>% as_tibble(rownames = "gene_id"), file.path(results_dir, paste0(prefix, "_detectedGenesNormalizedCounts_min_10_reads_in_one_condition.tsv")))

# Set the reference to the contrast level 2 (baseline) given by the --c2 option
dds[[cond_col]] = relevel( dds[[cond_col]], contrast[[3]])

# run DESeq
dds <- DESeq(dds, parallel = (n_cpus > 1))

# get normalized counts
nc <- counts(dds, normalized=T)

### IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds, filterFun=ihw, contrast=contrast) %>%
  as_tibble(rownames = "gene_id") %>%
  left_join(ensg_to_genesymbol) %>%
  left_join(ensg_to_desc, by = c("gene_id" = gene_id_type) ) %>%
  rename(genes_description = GENENAME) %>%
  arrange(pvalue)
summary(resIHW)
sum(resIHW$padj < fdr_cutoff, na.rm=TRUE)


# Filter for adjusted p-value < fdr_cutoff
resIHWsig <- resIHW %>% filter(padj < fdr_cutoff)

# significant genes as DE gene FDR < fdr_cutoff & abs(logfoldchange) > fc_cutoff , all genes as background
resIHWsig_fc <- resIHWsig %>% filter(abs(log2FoldChange) > fc_cutoff)

# Stop here if we do not have any DE genes
if(nrow(resIHWsig_fc) < 0) {
  stop("NO significant DE genes found: check fc_cutoff and fdr_cutoff!")
}

###### Perform Biotype QC
if(!is.null(gtf_file)) {
  gtf = rtracklayer::import(gtf_file, feature.type="gene") %>%
    as_tibble() %>%
    mutate(gene_id = remove_ensg_version(gene_id))

  count_before = nrow(resIHW)
  resIHW = resIHW %>% left_join(select(gtf, gene_id, gene_type), by=c("gene_id"="gene_id"))
  stopifnot("Number of genes should be the same after adding biotypes"= count_before == nrow(resIHW))
  resIHWsig = resIHWsig %>% left_join(select(gtf, gene_id, gene_type), by=c("gene_id"="gene_id"))

  biotype_counts = resIHWsig %>% group_by(gene_type) %>% count()

  p = biotype_counts %>%
    ggplot(aes(x=gene_type, y=n)) + geom_bar(stat='identity') + theme_bw() + coord_flip()
  ggsave(file.path(results_dir, paste0(prefix, "_biotype_counts.png")), p)
  write_tsv(biotype_counts, file.path(results_dir, paste0(prefix, "_biotype_counts.tsv")))

}

#### result list
de_res_list <- list(IHWallGenes = resIHW, IHWsigGenes = resIHWsig, IHWsigFCgenes = resIHWsig_fc)

#### write results to TSV and XLSX files
lapply(names(de_res_list), function(res) {
  fc_suffix <- ifelse(res == "IHWsigFCgenes", paste0("_", 2^fc_cutoff, "_fold"), "")
  write_tsv(de_res_list[[res]], file.path(results_dir, paste0(prefix, "_", res, fc_suffix, ".tsv")))
  write_xlsx(de_res_list[[res]], file.path(results_dir, paste0(prefix, "_" , res, fc_suffix, ".xlsx")))
})

###### Run TOPGO analysis
de_symbols <- resIHWsig_fc$gene_id
bg_symbols <- rownames(dds)[rowSums(counts(dds)) > 0]

lapply(c("BP", "MF"), function(ontology) {
  topgoDE <- topGOtable(de_symbols, bg_symbols,
                        ontology = ontology,
                        mapping = "org.Hs.eg.db",
                        geneID = gene_id_type)
  write_tsv(topgoDE, file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".tsv")))
  write_xlsx(topgoDE %>% select(-genes), file.path(results_dir, paste0(prefix, "_topGO_IHWsig_", ontology, ".xlsx")))
})


##### Pathway enrichment analysis
hgnc_to_entrez = AnnotationDbi::select(org.Hs.eg.db, resIHW %>% pull("gene_name") %>% unique(), keytype="SYMBOL", columns=c("ENTREZID"))

# full list with ENTREZIDs added
resIHW_entrez = resIHW %>%  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))
universe = resIHW_entrez %>% pull("ENTREZID") %>% unique()

# list of significant genes with ENTREZIDs added
resIHWsig_fc_entrez <- resIHWsig_fc %>%  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))
de_foldchanges <- resIHWsig_fc_entrez$log2FoldChange
names(de_foldchanges) <- resIHWsig_fc_entrez$ENTREZID

## ORA
ora_tests = list(
  "KEGG" = function(genes, universe) {
    enrichKEGG(
      gene         = genes,
      universe     = universe,
      organism     = 'hsa',
      pvalueCutoff = 0.05
    )
  },
  "Reactome" = function(genes, universe) {
    enrichPathway(
      gene = genes,
      organism = "human",
      universe = universe,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
  },
  "WikiPathway" = function(genes, universe) {
    enrichWP(
      gene = genes,
      universe     = universe,
      organism     = 'Homo sapiens',
      pvalueCutoff = 0.05
    )
  },
  "GO_BP" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  },
  "GO_MF" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = "MF",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  }
)

lapply(names(ora_tests), function(ora_name) {
  message(paste0("Performing ", ora_name, "ORA-test..."))
  test_fun = ora_tests[[ora_name]]
  ora_res = test_fun(resIHWsig_fc_entrez$ENTREZID, universe)
  ora_res = setReadable(ora_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  res_tab = as_tibble(ora_res@result)
  write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_ORA_", ora_name, ".tsv")))
  if (min(res_tab$p.adjust) < 0.05) {
    p = dotplot(ora_res, showCategory=40)
    ggsave(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_dotplot.png")), plot = p, width = 15, height = 10)

    p <- cnetplot(ora_res,
                  categorySize="pvalue",
                  showCategory = 5,
                  foldChange=de_foldchanges,
                  vertex.label.font=6)
    ggsave(file.path(results_dir, paste0(prefix, "_ORA_", ora_name, "_cnetplot.png")), plot = p, width = 15, height = 12)
  } else {
    message(paste0("Warning: No significant enrichment in ", ora_name, " ORA analysis. "))
  }
})

## GSEA
if(!skip_gsea) {
  # for GSEA use genes ranked by test statistic
  res_ihw_ranked = resIHW_entrez %>%
    arrange(-stat) %>%
    select(ENTREZID, stat) %>%
    na.omit() %>%
    distinct(ENTREZID, .keep_all=TRUE)
  ranked_gene_list = res_ihw_ranked$stat
  names(ranked_gene_list) = res_ihw_ranked$ENTREZID

  gsea_tests = list(
    "KEGG"=function(ranked_gene_list) {
      gseKEGG(geneList = ranked_gene_list, organism = "hsa", pvalueCutoff = 1)
    },
    "Reactome"=function(ranked_gene_list) {
      gsePathway(geneList = ranked_gene_list, organism = "human", pvalueCutoff = 1)
    },
    "WikiPathway"=function(ranked_gene_list) {
      gseWP(geneList = ranked_gene_list, organism = "Homo sapiens", pvalueCutoff = 1)
    },
    "GO_BP"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
               keyType = "ENTREZID",
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               minGSSize = 10)
    },
    "GO_MF"=function(ranked_gene_list) {
      gseGO(geneList=ranked_gene_list,
               keyType = "ENTREZID",
               OrgDb = org.Hs.eg.db,
               ont = "MF",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               minGSSize = 10)
    }
  )

  lapply(names(gsea_tests), function(gsea_name) {
    message(paste0("Performing ", gsea_name, " GSEA-test..."))
    test_fun = gsea_tests[[gsea_name]]
    gsea_res = test_fun(ranked_gene_list)
    gsea_res = setReadable(gsea_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    res_tab = gsea_res@result %>% as_tibble()

    write_tsv(res_tab, file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, ".tsv")))
    if (min(res_tab$p.adjust) < 0.05) {
      p = dotplot(gsea_res, showCategory=40)
      ggsave(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_dotplot.png")), plot = p, width = 15, height = 10)

      p <- cnetplot(gsea_res,
                    categorySize="pvalue",
                    showCategory = 5,
                    foldChange=de_foldchanges,
                    vertex.label.font=6)
      ggsave(file.path(results_dir, paste0(prefix, "_GSEA_", gsea_name, "_cnetplot.png")), plot = p, width = 15, height = 12)
    } else {
      message(paste0("Warning: No significant enrichment in ", gsea_name, " GSEA analysis. "))
    }
  })
}


########### PCA plot
vsd <- vst(dds, blind=FALSE)


p <- plotPCA(vsd, intgroup=c(cond_col)) +
  geom_point() +
  geom_text(vjust = 0,hjust = 0.2, nudge_x = -1, nudge_y = 0.5, aes(label = name)) +
  ggtitle(paste0("PCA: ", plot_title)) +
  scale_color_brewer(type="qual", palette="Set1") +
  theme_bw()
ggsave(file.path(results_dir, paste0(prefix, "_PCA.png")), plot = p)


########## Volcano plot
png(filename =  file.path(results_dir, paste0(prefix, "_volcano.png")), width = 650, height = 650, units = "px")
EnhancedVolcano(resIHW,
                lab = resIHW$gene_name,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 1e-6,
                FCcutoff = fc_cutoff,
                subtitle = "",
                legendPosition = "right",
                title = plot_title)
dev.off()

png(filename =  file.path(results_dir, paste0(prefix, "_volcano_padj.png")), width = 650, height = 650, units = "px")
EnhancedVolcano(resIHW,
                lab = resIHW$gene_name,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = fdr_cutoff,
                FCcutoff = fc_cutoff,
                subtitle = "",
                legendPosition = "right",
                title = plot_title)
dev.off()
