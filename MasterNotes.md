# Workflow/pipeline details for the albopictus autogeny RNA-seq analysis

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found in the following google sheets: ([link](<>))

### Data Accession
Data was generated from cateprillars reared at the University of North Carolina. Caterpillars were dissected and RNA was extracted by Sam Sturiale at Georgetown University.

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/trim.sh))

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/fastqc.sh))
