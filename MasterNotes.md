# Workflow/pipeline details for the Pieris rapae RNA-seq analysis

## Accessing Sequencing Data 
Maryland Genomics Facilty uses Aspera to transfer read files. So, I first had to download Aspera to the MAC desktop, then use Aspera to pull down the sequencing files using the following command then enter the password provided by UMD:

$HOME/insertfilepathhere/Applications/Aspera\ Connect.app/Contents/Resources/ascp -l 1024M -k 1 -QT insertlinkfromUMDhere $HOME/<filepath>

This resulted in three files: a tar filed containing the trimmed and untrimmed version of the reads for each sample (forward and reverse reads), an md5sum file, and a text summary file. 

The trimming done at UMD involved the following parameters: “We use Trimmomatic and trim for both adaptors and quality. We run Trimmomatic with parameters: simple clip threshold=7, seed mismatches=2, palindrome threshold=40, minimum sequence length=30, and training quality=20.” Also, from Lisa Sadzewicz: “Trimmed files included trimmed index/adaptor and quality trimming.” I ended up using the trimmed version but also doing my own trimming on top of the trimming done at UMD.

To check the md5sum file, I used the following commands and made sure the two outputs matched to ensure the tar file downloaded fully:
$md5sum insertnameoftarfilehere
$head insertnameofmd5sumfilehere 

Next, I needed to extract the tar file (do this on HPC so move tar file to google bucket then onto HPC. Use the upload.sh script for this). Then I used the extract.sh script to extract the tar file. Finally, I added the sample name to each of the file names. 


## Upstream 
Information about sequence reads (raw read count, read alignment rate, etc.) for the data can be found in the following google sheets: [link](https://docs.google.com/spreadsheets/d/1QJW4FL8r60wM-CdHDC86XxDxcfvSClsfac3ZJPZFG60/edit#gid=0)

### Data Accession
Data was generated from cateprillars reared at the University of North Carolina. Caterpillars were dissected and RNA was extracted by Sam Sturiale at Georgetown University.

The raw reads are available in NCBI’s short read archive (SRA) under accession number PRJNA1258715

### Preprocessing and Quality Control
Trimmomatic (version 0.39) was used to trim sequence reads based on quality ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/trim.sh))

FastQC (v0.11.9) was used for quality control visualization ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/fastqc.sh))

### Mapping
We will be mapping with STAR (v2.7.1a).

Mapping was done using the latest genome assembly available on NCBI for Pieris rapae (GCF_905147795.1 aka ilPieRapa1.1)

STAR (v2.7.1a) was used for indexing the genome ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/STAR_genomeIndex.sh))

Reads were mapped in a two pass method. The first pass followed typical method with splice junctions from annotations ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/STAR_mapping.sh)). The second pass is similar except that it additionally uses the output splice junctions info from the first pass (these would be novel splice junctions) to facilitate mapping ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/STAR_mapping_twopass.sh)). 

Output sam files were converted to bam (after which the sam files were deleted) and then bam files were indexed ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/sam2bam.sh)).

### Read Counting
HTSeq (v0.13.5) was used to counts reads mapped to genes for downstream analyses.
First used HTseq to check the strandedness. Illumina report from UMD specifies that they ran Strand-Specific RNA-sequencing, so the reads must be either forward or reverse stranded. To check which, I ran one file using forward, reverse, and no strandedness then chose the one with the fewest "no feature" outputs ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/htseq_stranded_test.sh)). For a test sample that began with 36,320,413 aligned reads, the forward run produced 36,009,968 no features while reverse run produced only 3,266,244 no features, so reads appear to be reverse. For more information about strandedness and running HTseq, see this paper: [link](https://academic.oup.com/bfg/article/19/5-6/339/5837822?login=false).

Then I used HTSeq to produce read count files for all samples ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/htseq.sh)).

## Downstream 

### DESeq2 (v1.44.0) package in R (v4.4.1) was used to identify differentially expressed genes ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/Prapae_DESeq_2factordesign.R))

### clusterProfiler (v4.12.6) packagein R was used to perform overrepresentation tests using both gene ontology (GO) terms and KEGG pathways ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/Prapae_Enrichment_analysis.R))

### the hclust function in R was used to perform utilized hierarchical clustering to examine the different patterns of expression across samples in the list of genes that were DE in response to the interaction of developmental and testing temperature ([script](https://github.com/samstur/P_rapae_RNASeq/blob/main/Prapae_IntGene_Clustering.R))

### This [script](https://github.com/samstur/P_rapae_RNASeq/blob/main/Prapae_Enrichment_Plot.R) was used to plot the enrichment results

### This [script](https://github.com/samstur/P_rapae_RNASeq/blob/main/Prapae_log2fc_rearing_v_testing_plot.R) was used to plot the figure showing the log2FC values for genes significantly affected by developmental temperature and short-term acclimation temperature.

### This [script](https://github.com/samstur/P_rapae_RNASeq/blob/main/Prapae_Mortality_DevTime_Plot.R) was used to analyze and plot data on mortality and developmental time for this experiment.

### This [script](https://github.com/samstur/P_rapae_RNASeq/blob/main/climate_UNC_fig.R) was used to download and process data on the climate of Chapel Hill, NC.

