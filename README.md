R.-ChIP-seq.histone
===================

Analysis of histone marks, and their differential presence in the genome

`Analysis.rs0.Rmd` - Using BED files of differential peaks defined by MACS2. Reads overlapping low complexity regions were removed by RepeatSoaker at 0 threshold. Peaks are annotated, checked for overlap, looked up where they are located, and tested for gene ontology and KEGG enrichment.

`Analysis.Rmd` - Using raw matrix counts and DESeq2 to compute the differences. Experimenting.

`Makefile` - NGS pipeline for processing BAM files.