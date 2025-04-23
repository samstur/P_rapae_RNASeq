# This script runs DESeq (specifically starting with htseq count files)
# a lot of this info is from https://github.com/hbctraining/DGE_workshop_salmon_online/blob/master/lessons/06_DGE_visualizing_results.md
########################## Installations ###########################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
BiocManager::install("apeglm")
BiocManager::install("GOplot")
BiocManager::install("mygene")
install.packages("RColorBrewer")
install.packages("ggolot2")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("viridis")
install.packages("vsn")
install.packages("ggthemes")
install.packages("VennDiagram")
BiocManager::install("genefilter")
BiocManager::install("ggrepel")
BiocManager::install("vsn")
BiocManager::install("ashr")
########################## Load Libraries ###########################
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(GOplot) 
library("apeglm")
library(mygene)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(vsn)
library(VennDiagram)
library(genefilter)
library(ggrepel)
library(ashr)
library(dplyr)
library(tidyr)
library(tibble)
library(mygene)
library(tidyr)
library(tidyverse)

########################## Input HTSeq data files ###########################
#Choose directory with htseq-count data
directory<-("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/HTSeq_output")

#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- list.files(directory)
head(sampleFiles)
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
head(sampleNames)
sampleNames <- sub(".*_L004_","",sampleNames) #this is removing the ending of the files to better represent the sample names
head(sampleNames)
sampleConditions <- substr(sampleNames, 1, 3) # this takes the first three letters of each sampleName to give the treatment name without the rep number
head(sampleConditions)
sampleConditions
sampleRearingTemps <- c("16","16","16","16","16","16","16","16",
                        "30","30","30","30","30","30","30","30")
sampleTestingTemps <- c("30","30","30","30","16","16","16","16",
                        "16","16","16","16","30","30","30","30")

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions,
                          RearingTemp = sampleRearingTemps,
                          TestingTemp = sampleTestingTemps) 

sampleTable$condition_new <- paste(sampleTable$RearingTemp,sampleTable$TestingTemp,sep="-")

sampleTable$condition <- factor(sampleTable$condition)
sampleTable$RearingTemp <- factor(sampleTable$RearingTemp)
sampleTable$TestingTemp <- factor(sampleTable$TestingTemp)
sampleTable$condition_new <- factor(sampleTable$condition_new)
View(sampleTable)

########################## Make the DESeq dataset from this HTSeq count data (this time with interaction term) ###########################

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ RearingTemp + TestingTemp + RearingTemp:TestingTemp)
dds

########################## Pre-filtering ###########################
#DESeq recommends a pre-filtering step to reduce memory size and increase speed. 
#They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

########################## Look at normalized data ###########################
# this isn't used downstream in the DE analysis because the DESeq() function does the normalization automatically behind the scenes
# instead we can use these normalized counts for downstream visualization if we want
dds_counts <- estimateSizeFactors(dds)
sizeFactors(dds_counts) # tells us the normalization factors for each sample
#counts(dds, normalized=TRUE) is providing counts scaled by size or normalization factors. 
dds_counts <- counts(dds_counts, normalized = TRUE) #save this to an excel file to look at normalized counts for visualization
head(dds_counts)
# Transform normalized counts using the vsd transformation 
cds <- estimateSizeFactors(dds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds, blind=TRUE)
theme_set(theme_bw())
meanSdPlot(assay(vsd))
plotDispEsts(cds)

data <- DESeq(dds)

## Extract raw counts that are normalized by size factor for later clustering analysis 
raw_counts <- counts(data, normalized = TRUE)
# first normalize these counts using vst function
norm_counts <- vst(data, blind = TRUE)%>% assay()
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")
raw_counts %>%
  as_tibble(rownames = "gene") %>%
  readr::write_csv("counts_raw.csv")
norm_counts %>%
  as_tibble(rownames = "gene") %>%
  readr::write_csv("counts_transformed.csv")

resultsNames(data)
res_rearingtemp <- results(data,name="RearingTemp_30_vs_16", alpha = 0.05)
res_testingtemp <- results(data,name="TestingTemp_30_vs_16", alpha = 0.05)
res_interaction <- results(data,name="RearingTemp30.TestingTemp30", alpha = 0.05)

res_rearingtemp.shrink <- lfcShrink(data,  
                                    type = "ashr",
                                    res=res_rearingtemp)
res_testingtemp.shrink <- lfcShrink(data, 
                                    type = "ashr",
                                    res=res_testingtemp)
res_interaction.shrink <- lfcShrink(data,  
                                    type = "ashr",
                                    res=res_interaction)


setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")
write.csv(file="./2factordesign_interaction_allgenes.csv",res_interaction.shrink, row.names = T)

res_rearingtemp.shrink_tbl <- res_rearingtemp.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_testingtemp.shrink_tbl <- res_testingtemp.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_interaction.shrink_tbl <- res_interaction.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 1
# Subset the significant results
res_rearingtemp.shrink_sig <- dplyr::filter(res_rearingtemp.shrink_tbl, 
                                            padj < padj.cutoff,
                                            abs(log2FoldChange) > lfc.cutoff)
res_testingtemp.shrink_sig <- dplyr::filter(res_testingtemp.shrink_tbl, 
                                            padj < padj.cutoff,
                                            abs(log2FoldChange) > lfc.cutoff)
res_interaction.shrink_sig <- dplyr::filter(res_interaction.shrink_tbl, 
                                            padj < padj.cutoff)

summary(res_rearingtemp.shrink_sig) #232 genes
summary(res_testingtemp.shrink_sig) #200 genes
summary(res_interaction.shrink_sig) # 51 genes if there is just a adj-p value cutoff and no lfc cutoff

# add gene names from NCBI to each of these data frames
library(mygene)
genes.rearingtemp <- queryMany(res_rearingtemp.shrink_sig$gene, scopes="symbol", fields=c("name"))
colnames(genes.rearingtemp)<-c("gene","id","score","gene name")
head(genes.rearingtemp)
genes.rearingtemp <- subset(genes.rearingtemp, select = c("gene", "gene name"))
res_rearingtemp.shrink_sig.gene <- as.data.frame(merge(genes.rearingtemp,res_rearingtemp.shrink_sig,by="gene"))

genes.testingtemp <- queryMany(res_testingtemp.shrink_sig$gene, scopes="symbol", fields=c("name"))
colnames(genes.testingtemp)<-c("gene","id","score","gene name")
head(genes.testingtemp)
genes.testingtemp <- subset(genes.testingtemp, select = c("gene", "gene name"))
res_testingtemp.shrink_sig.gene <- as.data.frame(merge(genes.testingtemp,res_testingtemp.shrink_sig,by="gene"))

genes.interaction <- queryMany(res_interaction.shrink_sig$gene, scopes="symbol", fields=c("name"))
colnames(genes.interaction)<-c("gene","id","score","gene name")
head(genes.interaction)
genes.interaction <- subset(genes.interaction, select = c("gene", "gene name"))
res_interaction.shrink_sig.gene <- as.data.frame(merge(genes.interaction,res_interaction.shrink_sig,by="gene"))

# save each list to a csv 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")
write.csv(res_rearingtemp.shrink_sig.gene, 
          file="./2factordesign_rearingtemp_30v16_sig.gene.csv", row.names = F)
write.csv(res_testingtemp.shrink_sig.gene, 
          file="./2factordesign_testingtemp_30v16_sig.gene.csv", row.names = F)
write.csv(res_interaction.shrink_sig.gene, 
          file="./2factordesign_interaction_sig.gene.csv", row.names = F)



#################################### Now again without the interaction ###########################################

dds_noint <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        design = ~ RearingTemp + TestingTemp)
dds_noint


########################## Pre-filtering ###########################
#DESeq recommends a pre-filtering step to reduce memory size and increase speed. 
#They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds_noint)) >= 10
dds_noint <- dds_noint[keep,]

data_noint <- DESeq(dds_noint)

resultsNames(data_noint)
res_rearingtemp_noint <- results(data_noint,name="RearingTemp_30_vs_16", alpha = 0.05)

res_testingtemp_noint <- results(data_noint,name="TestingTemp_30_vs_16", alpha = 0.05)

res_rearingtemp.shrink_noint <- lfcShrink(data_noint, 
                                          type = "ashr",
                                          res=res_rearingtemp_noint)
res_testingtemp.shrink_noint <- lfcShrink(data_noint,
                                          type = "ashr",
                                          res=res_testingtemp_noint)


## Convert to tibbles then subset the gene lists by ones that are significant 
res_rearingtemp.shrink_tbl_noint <- res_rearingtemp.shrink_noint %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_testingtemp.shrink_tbl_noint <- res_testingtemp.shrink_noint %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

## save this full list of genes with p-values to use later for KEGG enrichment
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")
write.csv(res_rearingtemp.shrink_noint, 
          file="./2factordesign_rearingtemp_30v16_noint_allgenes.csv", row.names = T)
write.csv(res_testingtemp.shrink_noint, 
          file="./2factordesign_testingtemp_30v16_noint_allgenes.csv", row.names = T)


# Now set thresholds for significantly differentially expressed genes 
padj.cutoff <- 0.05
lfc.cutoff <- 1
# Subset the significant results
res_rearingtemp.shrink_sig_noint <- dplyr::filter(res_rearingtemp.shrink_tbl_noint, 
                                                  padj < padj.cutoff,
                                                  abs(log2FoldChange) > lfc.cutoff)
res_testingtemp.shrink_sig_noint <- dplyr::filter(res_testingtemp.shrink_tbl_noint, 
                                                  padj < padj.cutoff,
                                                  abs(log2FoldChange) > lfc.cutoff)

summary(res_rearingtemp.shrink_sig_noint) #230 genes ## same number as if you run it with the reverse design
summary(res_testingtemp.shrink_sig_noint) #192 genes ## same number as if you run it with the reverse design

## remove the genes from the interaction list
res_rearingtemp.shrink_sig_noint <- subset(res_rearingtemp.shrink_sig_noint,
                                           !(res_rearingtemp.shrink_sig_noint$gene %in% res_interaction.shrink_sig.gene$gene))
summary(res_rearingtemp.shrink_sig_noint) # now 221 without the interaction genes

res_testingtemp.shrink_sig_noint <- subset(res_testingtemp.shrink_sig_noint,
                                           !(res_testingtemp.shrink_sig_noint$gene %in% res_interaction.shrink_sig.gene$gene))
summary(res_testingtemp.shrink_sig_noint) # now 185 without the interaction genes

# add gene names from NCBI to each of these data frames
genes.rearingtemp_noint <- queryMany(res_rearingtemp.shrink_sig_noint$gene, scopes="symbol", fields=c("name"))
colnames(genes.rearingtemp_noint)<-c("gene","id","score","gene name")
head(genes.rearingtemp_noint)
genes.rearingtemp_noint <- subset(genes.rearingtemp_noint, select = c("gene", "gene name"))
res_rearingtemp.shrink_sig.gene_noint <- as.data.frame(merge(genes.rearingtemp_noint,res_rearingtemp.shrink_sig_noint,by="gene"))
View(res_rearingtemp.shrink_sig.gene_noint)

genes.testingtemp_noint <- queryMany(res_testingtemp.shrink_sig_noint$gene, scopes="symbol", fields=c("name"))
colnames(genes.testingtemp_noint)<-c("gene","id","score","gene name")
head(genes.testingtemp_noint)
genes.testingtemp_noint <- subset(genes.testingtemp_noint, select = c("gene", "gene name"))
res_testingtemp.shrink_sig.gene_noint <- as.data.frame(merge(genes.testingtemp_noint,res_testingtemp.shrink_sig_noint,by="gene"))
View(res_testingtemp.shrink_sig.gene_noint)

summary(res_rearingtemp.shrink_sig.gene_noint) #221 genes
summary(res_testingtemp.shrink_sig.gene_noint) #185 genes 


## save each to csv
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")
write.csv(res_rearingtemp.shrink_sig.gene_noint, 
          file="./2factordesign_rearingtemp_30v16_sig.gene_noint.csv", row.names = F)
write.csv(res_testingtemp.shrink_sig.gene_noint, 
          file="./2factordesign_testingtemp_30v16_sig.gene_noint.csv", row.names = F)


### save both sets of DESeq results to RData file 
# Save my workspace to complete_image.RData in th,e
#  data folder of my working directory
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")
#save.image(file = "DESeq2_2factoranalysis.RData")

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles")

load("DESeq2_2factoranalysis.RData") 

sessionInfo()
