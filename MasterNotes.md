# Workflow/pipeline details for the Pieris rapae RNA-seq analysis

## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found in the following google sheets: ([link](<>))

### Data Accession
Data was generated from cateprillars reared at the University of North Carolina. Caterpillars were dissected and RNA was extracted by Sam Sturiale at Georgetown University.

The raw reads are available in NCBI’s short read archive (SRA) under accession number XXXXYYYY

### Preprocessing and Quality Control
Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/trim.sh))

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/fastqc.sh))

### Mapping
We will be mapping with STAR (v2.7.1a).

Mapping was done using the latest genome assembly available on NCBI for Pieris rapae (GCF_905147795.1 aka ilPieRapa1.1)

STAR (v2.7.1a) was used for indexing the genome ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/STAR_genomeIndex.sh))

Reads were mapped in a two pass method. The first pass followed typical method with splice junctions from annotations ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/STAR_mapping.sh)). The second pass is similar except that it additionally uses the output splice junctions info from the first pass (these would be novel splice junctions) to facilitate mapping ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/STAR_mapping_twopass.sh)). 

Output sam files were converted to bam (after which the sam files were deleted) and then bam files were indexed ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/sam2bam.sh))

## Read counting
HTSeq (v0.13.5) was used to counts reads mapped to genes for downstream analyses.
First used HTseq to check the strandedness. Ran one file using forward, reverse, and no strandedness then chose the one with the fewest "no feature" outputs ([script]https://github.com/samstur/P_rapae_RNASeq/blob/main/htseq_stranded_test.sh). 
Then used HTSeq to produce read count files for all samples ([script]https://github.com/samstur/P_rapae_RNASeq/blob/main/htseq.sh).

