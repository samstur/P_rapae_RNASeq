

library(ggplot2)
library(tidyr)
library(dplyr)
library("ggrepel")
library(viridis)

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
rearingtemptrt<- read.csv("2factordesign_rearingtemp_30v16_noint_allgenes.csv")
testingtemptrt<-read.csv("2factordesign_testingtemp_30v16_noint_allgenes.csv")

## now get rid of ones that are in the interaction list
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
interaction <- read.csv("2factordesign_interaction_sig.gene.csv")

length(rearingtemptrt$gene) #12298
rearingtemptrt_genes <- rearingtemptrt$gene[!rearingtemptrt$gene %in% interaction$gene]
length(rearingtemptrt_genes) # 12247
rearingtemptrt <- subset(rearingtemptrt, (rearingtemptrt$gene %in% rearingtemptrt_genes))
length(rearingtemptrt$gene) #12247

testingtemptrt_genes <- testingtemptrt$gene[!testingtemptrt$gene %in% interaction$gene]
length(testingtemptrt_genes)
testingtemptrt <- subset(testingtemptrt, (testingtemptrt$gene %in% testingtemptrt_genes))

rearingtemptrt<-rearingtemptrt[colnames(rearingtemptrt) %in% colnames(testingtemptrt),]

##combine the lists of all genes for rearing and testing temperature
combineddata<-merge(rearingtemptrt,testingtemptrt,by='gene')
View(combineddata)
## remove genes with NA values for both padj columns
combineddata<-combineddata[!is.na(combineddata$padj.x) & !is.na(combineddata$padj.y),] ## columns ending in .x are rearingtemp results while columns ending in .y are testingtemp
combineddata$sig<-"N"
combineddata[combineddata$padj.x<0.05 & abs(combineddata$log2FoldChange.x)>1,]$sig <- "rearingtemp"
combineddata[combineddata$padj.y<0.05 & abs(combineddata$log2FoldChange.y)>1,]$sig <- "testingtemp"
combineddata[combineddata$padj.y<0.05 & abs(combineddata$log2FoldChange.y)>1 & combineddata$padj.x<0.05 & abs(combineddata$log2FoldChange.x)>1,]$sig <- "both"

max(combineddata$log2FoldChange.x) # 7.53
max(combineddata$log2FoldChange.y) # 5.59
min(combineddata$log2FoldChange.x) # -3.6
min(combineddata$log2FoldChange.y) # -2.5

g.trt.overlap_l2fc1<-ggplot(combineddata,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_hline(yintercept=0,linetype='dashed')+
  geom_vline(xintercept=0,linetype='dashed')+
  geom_point(data=combineddata[combineddata$sig=="N",],aes(x=log2FoldChange.x,y=log2FoldChange.y),
             colour='darkgrey',alpha=0.2)+
  geom_point(data=combineddata[!combineddata$sig=="N",],aes(x=log2FoldChange.x,y=log2FoldChange.y,
                                                            fill=sig,colour=sig,size=sig),shape=21)+
  theme_bw()+
  xlab("Developmental Temperature Log2FC")+
  ylab('Short-term Acclimation Temperature Log2FC')+
  coord_cartesian(xlim = c(-4, 8), ylim = c(-3, 6))+
  scale_x_continuous(breaks=seq(from=-10,to=10,by=1))+
  scale_y_continuous(breaks=seq(from=-10,to=10,by=1))+
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=11),
        legend.text = element_text(size=9),
        legend.title= element_text(size=9))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("DE for Both Factors (N=12)", "DE for Developmental Temp. (N=209)",
                              "DE for Short-term Acc. Temp (N=173)"))+
 # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
#                    labels = c("DE for Both Factors (N=12)", "DE for Developmental Temp. (N=209)",
#                               "DE for Short-term Acc. Temp. (N=173)"))+
  scale_colour_manual(values = c("#ab7702", "#3a7a9e", "#017556"),
                    labels = c("DE for Both Factors (N=12)", "DE for Developmental Temp. (N=209)",
                               "DE for Short-term Acc. Temp (N=173)"))+
  scale_size_manual(values = c(1.7,1.2,1.2))+
  theme(legend.title = element_blank(),
        legend.position = c(0.76, 0.15),
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.key.size = unit(0.42, "cm"),
        legend.key.width = unit(0.42,"cm"),
        legend.spacing.y = unit(0, "pt"),
        legend.box.background = element_rect(color="black", size=0.35))+
  guides(fill = guide_legend(override.aes = list(shape=21,size=2)),
         shape = "none",size="none")
g.trt.overlap_l2fc1

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("trt_rearingtemp_testingtemp_overlap_l2fc1_colorblindfriendly1.png",plot=g.trt.overlap_l2fc1,
       dpi=600,height=4.5,width=6)

