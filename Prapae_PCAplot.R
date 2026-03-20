
## load Rdata file from DESeq2 script
load("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_DESeq2_2factordesign_outputfiles/DESeq2_2factoranalysis.RData")

library(vsn)
library(DESeq2)
library(ggplot2)

########################## PCA Visualization ##########################
# Transform normalized counts using the vsd transformation 
cds <- estimateSizeFactors(dds)
cds <- estimateDispersions(cds)
vsd = varianceStabilizingTransformation(cds, blind=TRUE)
theme_set(theme_bw())
meanSdPlot(assay(vsd))
plotDispEsts(cds)
PCA_data <- plotPCA(vsd, intgroup = c("condition_new","RearingTemp","TestingTemp"), returnData = TRUE)
percentVar <- round(100 * attr(PCA_data, "percentVar"))
nudge <- position_nudge(y = 4)

pcaplot<-ggplot(PCA_data, aes(x = PC1, y = PC2, position_nudge(y=10),label=condition_new)) +
  geom_point(size =3, position = position_jitter(w=0.05, h=0.05), aes(shape=condition_new,color=condition_new)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  # geom_text_repel(size=5)+
  scale_color_manual(values = c("#1E88E5","#D81B60","#1E88E5","#D81B60"))+
  scale_shape_manual(values= c(17,17,19, 19))+
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title= element_text(size=12),
        legend.key=element_blank())+
  guides(colour = guide_legend(title="Treatment"),
         shape = guide_legend(title="Treatment"))
# By default plotPCA() uses the top 500 most variable genes. 
# You can change this by adding the ntop= argument and specifying how many of the genes you want the function to consider.
pcaplot
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("PCAplot.png",plot=pcaplot,dpi=600,units='in',height=4.5,width=6)
