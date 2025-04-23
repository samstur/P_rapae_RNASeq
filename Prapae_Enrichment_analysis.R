install.packages("BiocManager")
BiocManager::install("org.Dm.eg.db")
library("org.Dm.eg.db")
library("tidyverse")
library(mygene)
BiocManager::install("clusterProfiler")
library("clusterProfiler")
library(GO.db)
BiocManager::install("GOfuncR")
BiocManager::install("GOSim")
library("AnnotationDbi")
library(reshape2)
install.packages("protti")
library("protti")
library(stringr)
library("GOSim")

## ## Read in full DESeq2 output
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
rearingtemp_all<- read.csv("2factordesign_rearingtemp_30v16_noint_allgenes.csv")
length(rearingtemp_all$gene)
testingtemp_all<-read.csv("2factordesign_testingtemp_30v16_noint_allgenes.csv")
interaction_all <- read.csv("2factordesign_interaction_allgenes.csv")

View(rearingtemp_all)

### Read in just DE genes 
rearingtemp_sig <- read.csv("2factordesign_rearingtemp_30v16_sig.gene_noint.csv")
testingtemp_sig <- read.csv("2factordesign_testingtemp_30v16_sig.gene_noint.csv")
interaction_sig <- read.csv("2factordesign_interaction_sig.gene.csv")
length(testingtemp_sig$gene) #185
length(rearingtemp_sig$gene) #221
length(interaction_sig$gene) #51


## get rid of genes with NA values for padj (thrown out by DESeq2)
rearingtemp_all <- rearingtemp_all %>%
  filter(!is.na(padj))
testingtemp_all <- testingtemp_all %>%
  filter(!is.na(padj))
interaction_all <- interaction_all %>%
  filter(!is.na(padj))

## read in Drosophila best blastx protein hit for each prapae gene 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/")
DM <- read.csv("prapae_allgene_blastx_uniq.csv")
head(DM)
length(DM$gene) # 9358
length(unique(DM$gene)) # 9358
hist(DM$pident)
hist(DM$evalue)
hist(DM$bitscore)
mean()hist(DM$qscov)

# Set thresholds
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/#:~:text=The%20bit%2Dscore%20provides%20a,50%20is%20almost%20always%20significant.
#pident.cutoff <- 10
#qscov.cutoff <- 10
bit.cutoff <- 45
evalue.cutoff <- 0.001
# Subset the significant results
DM2 <- dplyr::filter(DM, bitscore > bit.cutoff, evalue < evalue.cutoff)
#DM2 <- dplyr::filter(DM, evalue < evalue.cutoff)
length(DM2$gene) # with only an evalue cutoff of 0.001, I retain 9348 out of 9358 blast hits
## bit cutoff of 45 gets you 9223
# with a bitscore cutoff of 50 and a p-value cutoff of 0.001, I retain 8906 blastx results. If I reduce the bitscore cutoff to 40, I keep 9329 genes
View(DM2)
head(DM2)


genes.DM2 <- queryMany(DM2$drosophila_protein, scopes="Refseq", fields=c("symbol"))
colnames(genes.DM2)<-c("drosophila_protein","id","score","drosophila_gene")
head(genes.DM2)
genes.DM2 <- subset(genes.DM2, select = c("drosophila_protein", "drosophila_gene"))
DM2_gene <- as.data.frame(merge(DM2,genes.DM2,by="drosophila_protein"))

DM2_gene <- DM2_gene[!duplicated(DM2_gene), ]
head(DM2_gene)
length(unique(DM2_gene$gene)) #9348 ## updated 9223 for bitscore cutoff of 45
length(unique(DM2_gene$drosophila_gene)) #7177. ## updated 7105
View(DM2_gene)

# https://stackoverflow.com/questions/65180859/select-function-for-annotations-from-annotation-db-is-not-working
columns(org.Dm.eg.db)
help(columns, package="AnnotationDbi")


zz <- AnnotationDbi::select(org.Dm.eg.db, keys=keys(org.Dm.eg.db, "SYMBOL"),
                           columns=c("GO","SYMBOL","ONTOLOGY"), keytype="SYMBOL")
colnames(zz)<-c("drosophila_gene","go_id","evidence","ontology")
head(zz)
View(zz)


## for molecular function
zz_MF <- subset(zz,zz$ontology=="MF")
head(zz_MF)
zz_MF <- zz_MF[,c(1,2,4)]
zz_MF <- zz_MF[!duplicated(zz_MF), ]
head(zz_MF)
View(zz_MF)
length(unique(zz_MF$drosophila_gene)) #12393

# for biological process
zz_BP <- subset(zz,zz$ontology=="BP")
head(zz_BP)
zz_BP <- zz_BP[,c(1,2,4)]
zz_BP <- zz_BP[!duplicated(zz_BP), ]
View(zz_BP)

## add GO terms to list of homolog blastx results
DM2_gene_GO_MF <- merge(DM2_gene, zz_MF, by= "drosophila_gene")
DM2_gene_GO_BP <- merge(DM2_gene, zz_BP, by= "drosophila_gene")
View(DM2_gene_GO_BP)

length(unique(DM2_gene_GO_MF$gene)) #7930 ## need to figure out why I've gone from 9348 prapae genes in DM2_gene to 7930 here
## updated 7840
## possibly not every Drosophila gene has a GO term. I think it's that
head(setdiff(unique(DM2_gene$gene),unique(DM2_gene_GO_MF$gene)))

# make a dataframe with GO terms as first column and the matching prapae gene as the second column
## I'm not sure how to tell whether I have direct or indirect GO term annotations
## the tutorial says that if you are starting with direction annotation, you
## should use buildGOmap to infer indirect annotation and generate a dataframe to plug into enricher()
## so basically I have the most specific GO term for each gene, but I also want the ancestor GO terms for each gene ?
# https://bioconductor.org/packages/release/data/annotation/manuals/org.Dm.eg.db/man/org.Dm.eg.db.pdf
# "org.Dm.egGO is an R object that provides mappings between entrez gene identifiers and the GO
# identifiers that they are directly associated with. This mapping and its reverse mapping do NOT
# associate the child terms from the GO ontology with the gene. Only the directly evidenced terms
# are represented here"
head(DM2_gene_GO_MF2)
## For molecular function
DM2_gene_GO_MF2 <- DM2_gene_GO_MF[,c(8,3)]
DM2_gene_GO_MF2_new <- buildGOmap(DM2_gene_GO_MF2)
## add names to the GO IDs
gonames <- go2term(DM2_gene_GO_MF2_new$GO)
head(gonames)
DM2_gene_GO_MF2_new$term = gonames$Term[match(DM2_gene_GO_MF2_new$GO, as.character(gonames$go_id))]
head(DM2_gene_GO_MF2_new)
DM2_gene_GO_MF2_names <- DM2_gene_GO_MF2_new[,c(1,3)]
head(DM2_gene_GO_MF2_names)

# for biological processes
DM2_gene_GO_BP2 <- DM2_gene_GO_BP[,c(8,3)]
DM2_gene_GO_BP2_new <- buildGOmap(DM2_gene_GO_BP2)
## add names to the GO IDs
gonames <- go2term(DM2_gene_GO_BP2_new$GO)
head(gonames)
DM2_gene_GO_BP2_new$term = gonames$Term[match(DM2_gene_GO_BP2_new$GO, as.character(gonames$go_id))]
head(DM2_gene_GO_BP2_new)
DM2_gene_GO_BP2_names <- DM2_gene_GO_BP2_new[,c(1,3)]
head(DM2_gene_GO_BP2_names)


length(unique(DM2_gene_GO_MF2_new$Gene)) # 7602 #7840
length(unique(DM2_gene_GO_BP2_new$Gene)) # 7772 #8015


## write a function that let's me run everything at once (plots and tables)
GOFULL <- function(inputgenelist,universegenelist, term2genelist, term2namelist){
  data_GO <- enricher(
    gene = inputgenelist$gene,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universegenelist$gene,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE = term2genelist,
    TERM2NAME = term2namelist
  )
  data_GO_out = as.matrix(data_GO@result)
  data_GO_out <- as.data.frame(data_GO_out)
  data_GO_out$p.adjust <- as.numeric(data_GO_out$p.adjust)
  data_GO_out <- subset(data_GO_out,data_GO_out$p.adjust < 0.05)
  df <- data.frame("deg_up_list"=NA,"deg_up_names"=NA,"deg_down_list"=NA,"deg_down_names"=NA,"deg_up_num"=NA,"deg_down_num"=NA)
  head(df)
  for (i in 1:nrow(data_GO_out)) {
    deg_list <- (strsplit(sub("geneID","", data_GO_out[i,]["geneID"]), split = "/")[[1]])
    if (length(intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange > 1)$gene))>0){
      deg_up_list <- intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange > 1)$gene)
      deg_up_names <- paste(subset(inputgenelist,inputgenelist$gene %in% deg_up_list)$gene.name,collapse=",")
      deg_up_num <- length(deg_up_list)
      deg_up_list <- paste(deg_up_list,collapse=",")
    } else {deg_up_list <- NA
    deg_up_num <- 0
    deg_up_names <- NA}
    if (length(intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange < 1)$gene))>0){
      deg_down_list <- intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange < 1)$gene)
      deg_down_names <- paste(subset(inputgenelist,inputgenelist$gene %in% deg_down_list)$gene.name,collapse=",")
      deg_down_num <- length(deg_down_list)
      deg_down_list <- paste(deg_down_list,collapse=",")
    } else {
      deg_down_list <- NA
      deg_down_num <- 0
      deg_down_names <- NA}
    newdata <- data.frame(deg_up_list,deg_up_names,deg_down_list,deg_down_names,deg_up_num,deg_down_num)
    rownames(newdata) <- NULL
    df<-rbind(df,newdata)
  }
  df<-df[!(row.names(df) == 1),] #remove NA fake row
  out_final<-cbind(data_GO_out,df)
  return(out_final)
}

##########################################################################
##########################################################################
##########################################################################
# same but starting with GO terms that I got off of FTP through NCBI for p rapae
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output/GO Tests")
goprapae <- read.csv("GCF_905147795.1_ilPieRapa1.1_gene_ontology.csv")
View(goprapae)

## separate by BP and MF
goprapae_BP <- subset(goprapae,goprapae$Aspect=="P")
goprapae_BP <- goprapae_BP[,c(5,3)]
colnames(goprapae_BP)<-c("go_id","gene")
goprapae_BP <- goprapae_BP[!duplicated(goprapae_BP), ]
length(unique(goprapae_BP$gene)) #7132
goprapae_BP_new <- buildGOmap(goprapae_BP)
head(goprapae_BP_new)
## add names to the GO IDs
gonames <- go2term(goprapae_BP_new$go_id)## changed this column name from GO to go_id
head(gonames)
goprapae_BP_new$term = gonames$Term[match(goprapae_BP_new$go_id, as.character(gonames$go_id))] ## changed this column name from GO to go_id
head(goprapae_BP_new)
goprapae_BP_names <- goprapae_BP_new[,c(1,3)]
head(goprapae_BP_names)


goprapae_MF <- subset(goprapae,goprapae$Aspect=="F")
goprapae_MF <- goprapae_MF[,c(5,3)]
colnames(goprapae_MF)<-c("go_id","gene")
goprapae_MF <- goprapae_MF[!duplicated(goprapae_MF), ]
View(goprapae_MF)
length(unique(goprapae_MF$gene)) #8004
goprapae_MF_new <- buildGOmap(goprapae_MF)
## add names to the GO IDs
gonames <- go2term(goprapae_MF_new$GO)
head(gonames)
goprapae_MF_new$term = gonames$Term[match(goprapae_MF_new$GO, as.character(gonames$go_id))]
head(goprapae_MF_new)
goprapae_MF_names <- goprapae_MF_new[,c(1,3)]
head(goprapae_MF_names)


##########################################################################
## now want to combine GO annotations from interproscan NCBI and from Drosophila homology results

combined_BP <- rbind(goprapae_BP_new,DM2_gene_GO_BP2_new)
combined_BP <- combined_BP[!duplicated(combined_BP), ]
View(combined_BP) # 452,292 entries
View(goprapae_BP_new) #198,274 entries
View(DM2_gene_GO_BP2_new) # 417,746 entries
length(unique(combined_BP$Gene)) #9102
length(unique(goprapae_BP_new$Gene)) #7132
length(unique(DM2_gene_GO_BP2_new$Gene)) #8015

combined_BP_names <- rbind(goprapae_BP_names,DM2_gene_GO_BP2_names)
combined_BP_names <- combined_BP_names[!duplicated(combined_BP_names), ]


combined_MF <- rbind(goprapae_MF_new,DM2_gene_GO_MF2_new)
combined_MF <- combined_MF[!duplicated(combined_MF), ]
View(combined_MF) # 106,807 entries
View(goprapae_MF_new) #85,350 entries
View(DM2_gene_GO_MF2_new) # 90,921 entries
length(unique(combined_MF$Gene)) #9365
length(unique(goprapae_MF_new$Gene)) #8004
length(unique(DM2_gene_GO_MF2_new$Gene)) #7804

combined_MF_names <- rbind(goprapae_MF_names,DM2_gene_GO_MF2_names)
combined_MF_names <- combined_MF_names[!duplicated(combined_MF_names), ]

# Save my workspace to complete_image.RData in th,e
#  data folder of my working directory
#setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output/GO Tests")
#save.image(file = "GOAnalysis_combined.RData")
# Load objects into my workspace
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output/GO Tests")
load(file = "GOAnalysis_combined.RData")

library(clusterProfiler)
library(stringr)

testing_BP_combined <- GOFULL(testingtemp_sig,testingtemp_all, combined_BP, combined_BP_names)
rearing_BP_combined <- GOFULL(rearingtemp_sig,rearingtemp_all, combined_BP, combined_BP_names)
testing_MF_combined <- GOFULL(testingtemp_sig,testingtemp_all, combined_MF, combined_MF_names)
rearing_MF_combined <- GOFULL(rearingtemp_sig,rearingtemp_all, combined_MF, combined_MF_names)

interaction_BP_combined <- GOFULL(interaction_sig,interaction_all, combined_BP, combined_BP_names)
interaction_MF_combined <- GOFULL(interaction_sig,interaction_all, combined_MF, combined_MF_names)
View(interaction_BP_combined)
View(testing_BP_combined)

#setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output/GO Tests")
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Enrichment_analysis_outputfiles")

write.csv(testing_BP_combined,"GO_testing_BP_combined.csv",row.names = F)
write.csv(rearing_BP_combined,"GO_rearing_BP_combined.csv",row.names = F)
write.csv(testing_MF_combined,"GO_testing_MF_combined.csv",row.names = F)
write.csv(rearing_MF_combined,"GO_rearing_MF_combined.csv",row.names = F)

write.csv(interaction_MF_combined,"GO_interaction_MF_combined.csv",row.names = F)
write.csv(interaction_BP_combined,"GO_interaction_BP_combined.csv",row.names = F)





##########################################################################
##########################################################################
##########################################################################
## same but do overrepresentation with KEGG instead of GO terms

## write a function that let's me run everything at once (plots and tables)
KEGGFULL <- function(inputgenelist,universegenelist){
  inputgenelist$ncbi_geneid <- sub("LOC","", inputgenelist$gene)
  universegenelist$ncbi_geneid <- sub("LOC","", universegenelist$gene)
  kk <- enrichKEGG(gene=inputgenelist$ncbi_geneid, 
                   universe=universegenelist$ncbi_geneid,
                   organism="prap", 
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = "BH",
                   keyType = "ncbi-geneid")
  data_KEGG_out = as.matrix(kk@result)
  data_KEGG_out <- as.data.frame(data_KEGG_out)
  data_KEGG_out$p.adjust <- as.numeric(data_KEGG_out$p.adjust)
  data_KEGG_out <- subset(data_KEGG_out,data_KEGG_out$p.adjust < 0.05)
  df <- data.frame("deg_up_list"=NA,"deg_up_names"=NA,"deg_down_list"=NA,"deg_down_names"=NA,"deg_up_num"=NA,"deg_down_num"=NA)
  for (i in 1:nrow(data_KEGG_out)) {
    deg_list <- (strsplit(sub("geneID","", data_KEGG_out[i,]["geneID"]), split = "/")[[1]])
    if (length(intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange > 1)$ncbi_geneid))>0){
      deg_up_list <- intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange > 1)$ncbi_geneid)
      deg_up_names <- paste(subset(inputgenelist,inputgenelist$ncbi_geneid %in% deg_up_list)$gene.name,collapse=",")
      deg_up_num <- length(deg_up_list)
      deg_up_list <- gsub(" ","",paste("LOC",deg_up_list,collapse=","))
    } else {deg_up_list <- NA
    deg_up_num <- 0
    deg_up_names <- NA}
    if (length(intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange < 1)$ncbi_geneid))>0){
      deg_down_list <- intersect(deg_list,subset(inputgenelist,inputgenelist$log2FoldChange < 1)$ncbi_geneid)
      deg_down_names <- paste(subset(inputgenelist,inputgenelist$ncbi_geneid %in% deg_down_list)$gene.name,collapse=",")
      deg_down_num <- length(deg_down_list)
      deg_down_list <- gsub(" ","",paste("LOC",deg_down_list,collapse=","))
    } else {
      deg_down_list <- NA
      deg_down_num <- 0
      deg_down_names <- NA}
    newdata <- data.frame(deg_up_list,deg_up_names,deg_down_list,deg_down_names,deg_up_num,deg_down_num)
    rownames(newdata) <- NULL
    df<-rbind(df,newdata)
  }
  df<-df[!(row.names(df) == 1),] #remove NA fake row
  out_final<-cbind(data_KEGG_out,df)
  return(out_final)
}


testingtempKEGG <- KEGGFULL(testingtemp_sig,testingtemp_all)
rearingtempKEGG <- KEGGFULL(rearingtemp_sig,rearingtemp_all)
interactionKEGG <- KEGGFULL(interaction_sig,interaction_all)


View(testingtempKEGG)
View(rearingtempKEGG)
View(interactionKEGG)


setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Enrichment_analysis_outputfiles")
write.csv(testingtempKEGG,"testingtempKEGG.csv",row.names = F)
write.csv(rearingtempKEGG,"rearingtempKEGG.csv",row.names = F)
write.csv(interactionKEGG,"interactionKEGG.csv",row.names = F)



