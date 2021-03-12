# R script to run for "Bioinformatic Service Unit" projects


## Usage
```
runDESeq2_ICBI.R

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
  --paired_grp=<paired_grp>     Column that conatins the name of the paired samples, when dealing with
                                paired data.
  --plot_title=<title>          Title shown above plots. Is built from contrast per default.
  --prefix=<prefix>             Results file prefix. Is built from contrasts per default.
  --fdr_cutoff=<fdr>            False discovery rate for GO analysis and volcano plots [default: 0.1]
  --fc_cutoff=<log2 fc cutoff>  Fold change (log2) cutoff for volcano plots [default: 1]
  --gtf_file=<gtf>              Path to the GTF file used for featurecounts. If specified, a Biotype QC
                                will be performed.
  --gene_id_type=<id_type>      Type of the identifier in the `gene_id` column compatible with AnnotationDbi [default: ENSEMBL] 
```

The R script takes the merged count table from the nf-core rnaseq pipeline as input (featurecounts.merged.counts.tsv, see example).

As second input file you need the sample table in csv format. See "sampleTableN.tsv" as example. We suggested to use the sample table that was used for running the nf-core rnaseq pipeline.

| group | replicate | pair_group | fastq_1 | fastq_2 | strandedness |
| ------| --------- | ---------- | ------- | ------- | ------------ |
| WT | 1 | | wt1_1.fq.gz | wt1_2.fq.gz | forward |
| WT | 2 | | wt2_1.fq.gz | wt2_2.fq.gz | forward |
| WT | 3 | | wt3_1.fq.gz | wt3_2.fq.gz | forward |
| KO | 1 | | ko1_1.fq.gz | ko1_2.fq.gz | forward |
| KO | 2 | | ko2_1.fq.gz | ko2_2.fq.gz | forward |
| KO | 3 | | ko3_1.fq.gz | ko3_2.fq.gz | forward |


## Note:
You need to run the 3.x version of the RNA-seq pipeline. The current release of
3.0 contains a bug, so use the dev version currently! 

## Methods summary

E-mail template to send out to "clients". Make sure to double check 
and maybe update version numbers. 

```
 * We processed raw FASTQ fils using the nf-core RNA-seq pipeline [doi:10.5281/zenodo.1400710, doi:10.1038/s41587-020-0439-x] version 3.1. 
 * In brief, FASTQ files were trimmed using trimgalore v0.6.6 and reads were aligned to the GRCh38 reference
   genome with GENCODE v33 annotations using STAR v2.7.6a [doi:10.1093/bioinformatics/bts635]. 
   Gene expression was quantified using Salmon v1.4.0 [doi:10.1038/nmeth.4197]
   using the aligned BAM files as input. 
 * We performed differential gene expression analysis with DESeq2 v1.30.0 [doi: 10.1186/s13059-014-0550-8].
   False-discovery-rates were calculated using IHW v1.18.0 [doi:10.1038/nmeth.3885]. 
   We performed gene ontology (GO)-term enrichment analysis using the topGO package v2.42.0.
   Volcano plots were visualized using the EnhancedVolcano package v1.8.0. 
```


## Output description

### Over-representation test vs Gene-Set-Enrichment Analysis
The over-representation test (ORA) takes the list of differentially expressed genes (as defined by the FDR and fold-change cutoff)
and checks if it contains more genes from a certain gene-set (e.g. GO-term, KEGG pathway) than would be expected
by random chance. 

Gene-set-enrichment-analysis (GSEA) rankes genes from most upregulated to most downregulated, and checks if genes from a certain gene-set
tend to be more at the top, or more at the bottom of this list. 

The ORA test results are straightforward to interpret and, in most cases, the preferred flavor of gene-set analysis. 
GSEA is superior when very few genes are differentially expressed, as it also takes into account more subtle, 
but coordinated gene expression changes that are statistically significant only when regarding the gene set as a whole, 
but not on the level of single genes. 

### Differential gene expression
 
 * `IHWallGenes`: Fold change and pvalues for all genes
 * `IHWsigFCgenes`: All differentially expressed genes that meet a fold change of 2 and a FDR of 0.1
 * `IHWsigGenes`: All differentially expressed genes that meet a FDR of 0.1. 
 * `PCA`: Principal component analysis plot. Ideally samples cluster by condition. 
 * `biotype_counts`: Lists which type of genomic features are most differentially expressed. Usually, we expect most of them being protein coding. 
 * `detectedGenesNormalizedCounts`: The gene expression matrix normalized by library size. Genes that are not expressed with at least 10 reads in any condition are filtered out. 
 * `detectedGenesRawCounts`: The raw gene expression matrix. Genes that are not expressed with at least 10 reads in any condition are filtered out. 
 * `go_enrich`: Plots and table with Gene-ontology enrichment results (using ORA implemented in clusterProfiler). 
 * `kegg_enrich`: Plots and table with KEGG pathway enrichment results (using ORA implemented in clusterProfiler). 
 * `reactome_enrich`: Plots and table with Reactome pathway enrichment results (using ORA implemented in clusterProfiler). 
 * `go_gsea`: Plots and table with Gene set 
 * `topGO`: tables with GO-term enrichment analysis by GO category (using the topGO method)
 * `volcano`: Volcano plot of log-fold change against unadjusted pvalues. 
 * `volcano_padj`: Volcano plot of log-fold change against FDRs. 
 * `wp`: Table with WikiPathway enrichment results (using the clusterProfiler method). 
