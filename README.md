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


### Note:
You need to run the **dev version** of the nf-core rnaseq pipeline. Was tested with Gregor Sturms fork (**grst/rnaseq -r e4e55f5**)


