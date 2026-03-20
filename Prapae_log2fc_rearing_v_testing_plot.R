

library(ggplot2)
library(tidyr)
library(dplyr)
library("ggrepel")
library(viridis)
library(ggvenn)


setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
rearingtemptrt<- read.csv("2factordesign_rearingtemp_30v16_noint_allgenes.csv")
testingtemptrt<-read.csv("2factordesign_testingtemp_30v16_noint_allgenes.csv")

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
max(combineddata$log2FoldChange.y) # 5.59 ## new is 6.31
min(combineddata$log2FoldChange.x) # -3.6
min(combineddata$log2FoldChange.y) # -2.5

### now add a column for combineddata which says whether each genes is in the sig int gene list
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
interaction <- read.csv("2factordesign_interaction_sig.gene.csv")
combineddata$intgene <- ""
combineddata <- combineddata %>%
  mutate(intgene = case_when(gene %in% interaction$gene ~ "Yes",
                             !gene %in% interaction$gene ~ "No"))

## one of the 51 interaction genes is not in the combineddata set because DESeq2 removed it as an outlier
## and set its p-adj and p-value to NA. So it's not in this plot
# the gene is LOC123689103 
setdiff(interaction$gene, subset(combineddata,combineddata$intgene=="Yes")$gene)

table(combineddata$intgene) ## so only 50 interaction genes 
table(combineddata$sig)
table(combineddata$sig,combineddata$intgene)

g.trt.overlap_l2fc1<-ggplot(combineddata,aes(x=log2FoldChange.x,y=log2FoldChange.y))+
  geom_hline(yintercept=0,linetype='dashed')+
  geom_vline(xintercept=0,linetype='dashed')+
  geom_point(data=combineddata[combineddata$sig=="N",],aes(x=log2FoldChange.x,y=log2FoldChange.y,shape=intgene),
             colour='darkgrey',alpha=0.2)+
  geom_point(data=combineddata[!combineddata$sig=="N",],aes(x=log2FoldChange.x,y=log2FoldChange.y,
                                                            fill=sig,colour=sig,size=sig,shape=intgene))+
  theme_bw()+
  xlab(bquote("Developmental Temperature Log"[2]*FC))+
  ylab(bquote("Short-term Acclimation Temperature Log"[2]*FC))+
  coord_cartesian(xlim = c(-4, 9), ylim = c(-3, 6.5))+
  scale_x_continuous(breaks=seq(from=-10,to=10,by=1))+
  scale_y_continuous(breaks=seq(from=-10,to=10,by=1))+
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=11),
        legend.text = element_text(size=9),
        legend.title= element_text(size=9))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("Sig. Effect of Both Factors", "Sig. Effect of Dev. Temp.",
                              "Sig. Effect of Acc. Temp."))+
 # scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
#                    labels = c("DE for Both Factors (N=12)", "DE for Developmental Temp. (N=209)",
#                               "DE for Short-term Acc. Temp."))+
  scale_colour_manual(values = c("#ab7702", "#3a7a9e", "#017556"),
                    labels = c("Sig. Effect of Both Factors", "Sig. Effect of Dev. Temp.",
                               "Sig. Effect of Acc. Temp."))+
  scale_shape_manual(values = c(21,24),
                     labels = c("No Significant Interaction",
                                "Significant Interaction"))+
  scale_size_manual(values = c(1.7,1.2,1.2))+
  theme(legend.title = element_blank(),
        legend.position = c(0.77, 0.76),
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.key.size = unit(0.43, "cm"),
        legend.key.width = unit(0.43,"cm"),
        legend.spacing.y = unit(0, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_rect(fill = "white", color = "black", size=0.37))+
  guides(fill = guide_legend(override.aes = list(size=2)),size="none")+
  labs(tag="A") 
g.trt.overlap_l2fc1

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("trt_rearingtemp_testingtemp_overlap_l2fc1_colorblindfriendly1_includesintgenes.png",plot=g.trt.overlap_l2fc1,
       dpi=600,height=4.5,width=6)

View(combineddata)

### pearson correlation coefficient is the same as vector analysis results 
cor(combineddata[combineddata$sig!="N",]$log2FoldChange.y, combineddata[combineddata$sig!="N",]$log2FoldChange.x, method = "pearson")
# r = -0.037
cor(combineddata$log2FoldChange.y, combineddata$log2FoldChange.x, method = "pearson")
# 0.099

# same as above but without int genes 
cor(combineddata[combineddata$sig!="N" & combineddata$intgene=="No",]$log2FoldChange.y, combineddata[combineddata$sig!="N" & combineddata$intgene=="No",]$log2FoldChange.x, method = "pearson")
# -0.061
cor(combineddata[combineddata$intgene=="No",]$log2FoldChange.y, combineddata[combineddata$intgene=="No",]$log2FoldChange.x, method = "pearson")
# 0.082

##################################################################################
## Now make a Venn diagram that shows overlap between the two main effects and the interaction genes

## need to put all three gene lists in a list together 

devsig <- combineddata[combineddata$padj.x<0.05 & abs(combineddata$log2FoldChange.x)>1,]$gene
length(devsig)

accsig <- combineddata[combineddata$padj.y<0.05 & abs(combineddata$log2FoldChange.y)>1,]$gene
length(accsig)

intsig <- interaction$gene
intsig

x <- list(A = devsig, B = accsig, C = intsig)
names(x) <- c("Dev. Temp.","Short-term
Acc. Temp.","Interaction")

venndiagram <- ggvenn(x, 
  fill_color = c("#56B4E9", "#009E73", "white"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F) +
  labs(tag="B") +
  theme(
    plot.margin = margin(t = 0, r = 0, b = 0, l = 1.2, unit = "cm")
  )
venndiagram

library(patchwork)
finalplot <- g.trt.overlap_l2fc1 + venndiagram

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("log2fc_venndiagram_includesintgenes_plot.png",plot=finalplot,
       dpi=600,height=5,width=10)
