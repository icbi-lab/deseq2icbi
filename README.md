# R script to run for "Bioinformatic Service Unit" projects


## Usage
```
The R script takes the merged count table from the nf-core rnaseq pipeline as input (merged_gene_counts.txt).
To tidy up the sample names in the count file, please run the *format_merged_gene_counts.sh* to do so.
As second input file you need the sample table in tsv format. See "sampleTableN.tsv" as example

Please edit the the parameters section at the begining of the script to your needs. (This needs to be implemented cleaner as cmd line args)

```

### Note:
 If you run the dev version of the nf-core rnaseq pipeline, please make sure that you use the option for generating the merged counts file (featurecounts.merged.counts.tsv ).


