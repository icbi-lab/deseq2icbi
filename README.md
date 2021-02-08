# R script to run for "Bioinformatic Service Unit" projects


## Usage
```
Usage:
  runDESeq2_ICBI_p.R <sample_sheet> <count_table> --result_dir=<res_dir> --c1=<c1> --c2=<c2> [options]
  runDESeq2_ICBI_p.R --help

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


