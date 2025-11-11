

install.packages("PieGlyph")
library("PieGlyph")
library("ggplot2")
library(stringr)
library("cowplot")
library(gridExtra)

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Enrichment_analysis_outputfiles")
GO_mf_testing <- read.csv("GO_testing_MF_combined.csv")
GO_mf_rearing <- read.csv("GO_rearing_MF_combined.csv")

GO_mf_testing$treatment <- "testing"
GO_mf_rearing$treatment <- "rearing"
dfboth_mf <- rbind(GO_mf_rearing,GO_mf_testing)


## first subset list to get rid of redundant parent terms
## got this list by running each list of rearing and testing terms in Revigo (http://revigo.irb.hr/Results?jobid=225187156)
## with small (0.5 setting, provided adjusted p-values and said lower is better, used the whole uniprot default list, and used SimRel)
## then I grouped them by similarity in genes manually 

newlist_mf <- c("unfolded protein binding","protein folding chaperone","heat shock protein binding","misfolded protein binding",
                "flavin adenine dinucleotide binding","oxidoreductase activity, acting on CH-OH group of donors",
                "hydrolase activity, hydrolyzing O-glycosyl compounds","hydrolase activity, acting on glycosyl bonds",
                "lipase activity","triglyceride lipase activity",
                "structural constituent of chitin-based cuticle",
                "carbon-nitrogen lyase activity","racemase and epimerase activity",
                "peptidase inhibitor activity")

dfboth_mf <- dfboth_mf[dfboth_mf$Description %in% newlist_mf,]


dfboth_mf_sort <- with(dfboth_mf,  dfboth_mf[order(Description) , ])

dfboth_mf_sort$Description
dfboth_mf_sort$linegroup <- c("lyase1","flavin2","HSP1","hydrolase1","hydrolase1",
                              "lipid2","HSP1","flavin2","peptidase","HSP2","HSP1",
                              "lyase1","cuticle1","lipid2","HSP2","HSP1")



mfplot<-ggplot(data=dfboth_mf_sort,aes(x=treatment,y=Description,group=linegroup))+
  geom_line(size=4,color="lightgrey")+
  geom_pie_glyph(slices = c('deg_up_num', 'deg_down_num'),radius=0.4)+
  geom_point(size=5,colour="white")+
  theme_bw()+
  xlab("")+
  ylab("")+
 # labs(tag="A")+
  ggtitle("(A) GO Molecular Function")+
  scale_x_discrete(limits=c("rearing","testing"),labels=c("Dev. 
Temp.","Short-term Acc. 
Temp."))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25),
                   limits = c("unfolded protein binding","protein folding chaperone","heat shock protein binding","misfolded protein binding",
                              "flavin adenine dinucleotide binding","oxidoreductase activity, acting on CH-OH group of donors",
                              "hydrolase activity, hydrolyzing O-glycosyl compounds","hydrolase activity, acting on glycosyl bonds",
                              "lipase activity","triglyceride lipase activity",
                              "structural constituent of chitin-based cuticle",
                              "carbon-nitrogen lyase activity","racemase and epimerase activity",
                              "peptidase inhibitor activity"))+
  scale_fill_manual(values=c("#EC0000","#0055BE"),name = "",labels=c("Upregulated","Downregulated"))+
  #scale_radius(breaks = c(3, 6, 9, 12), name="Number of DEG")+
  geom_label(aes(label=Count),fill=NA,label.size=NA,size=3)+
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title= element_text(size=10),
        plot.title = element_text(size = 11, face = "bold"),
        legend.margin=margin(t = 0, unit='cm'),
        legend.position = "none")
#        legend.spacing.y = unit(0.5, 'cm'))
mfplot

######### now for biological process
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Enrichment_analysis_outputfiles")
GO_bp_testing <- read.csv("GO_testing_bp_combined.csv")
GO_bp_rearing <- read.csv("GO_rearing_bp_combined.csv")

GO_bp_testing$treatment <- "testing"
GO_bp_rearing$treatment <- "rearing"
dfboth_bp <- rbind(GO_bp_rearing,GO_bp_testing)

dfboth_bp_sort <- with(dfboth_bp,  dfboth_bp[order(Description) , ])

## edited list to remove redundant parent terms 
finallist <- c("protein folding","protein refolding","regulation of translational initiation by eIF2 alpha phosphorylation","regulation of translation in response to stress",
               "xenobiotic metabolic process", "response to xenobiotic stimulus",
               "epidermal cell differentiation","cuticle development",
               "cellular ketone metabolic process","hormone metabolic process",
               "branched-chain amino acid metabolic process")


## first subset list to get rid of redundant parent terms
dfboth_bp_sort_new <- dfboth_bp_sort[dfboth_bp_sort$Description %in% finallist,]

dfboth_bp_sort_new$Description

linegroup_bp <- c("amino1","hormone2","epi2","epi2","hormone2","HSP1","HSP2","HSP1","HSP2","HSP1",
                  "HSP2","HSP1","xeno2","xeno2")
dfboth_bp_sort_new$linegroup_bp <- linegroup_bp


bpplot<-ggplot(data=dfboth_bp_sort_new,aes(x=treatment,y=Description,group=linegroup_bp))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25),
                   limits=c("protein folding","protein refolding","regulation of translational initiation by eIF2 alpha phosphorylation","regulation of translation in response to stress",
                            "xenobiotic metabolic process", "response to xenobiotic stimulus",
                            "epidermal cell differentiation","cuticle development",
                            "cellular ketone metabolic process","hormone metabolic process",
                            "branched-chain amino acid metabolic process"))+
  scale_fill_manual(values=c("#EC0000","#0055BE"),name = "",labels=c("Upregulated","Downregulated"))+
  geom_line(size=4,color="lightgrey")+
  geom_pie_glyph(slices = c('deg_up_num', 'deg_down_num'),radius=0.4)+
  geom_point(size=5,colour="white")+
  geom_label(aes(label=Count),fill=NA,label.size=NA,size=3)+
  theme_bw()+
  xlab("")+
  ylab("")+
  #labs(tag="B")+
  ggtitle("(B) GO Biological Process")+
  scale_x_discrete(limits=c("rearing","testing"),labels=c("Dev. 
Temp.","Short-term Acc. 
Temp."))+
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title= element_text(size=10),
        plot.title = element_text(size = 11, face = "bold"),
        #        legend.margin=margin(t = 0, unit='cm'),
        legend.position = "")
#        legend.spacing.y = unit(0.5, 'cm'))
bpplot

bpplot_legend<-ggplot(data=dfboth_bp_sort_new,aes(x=treatment,y=Description,group=linegroup_bp))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25),
                   limits=c("protein folding","protein refolding","regulation of translational initiation by eIF2 alpha phosphorylation","regulation of translation in response to stress",
                            "xenobiotic metabolic process", "response to xenobiotic stimulus",
                            "epidermal cell differentiation","cuticle development",
                            "cellular ketone metabolic process","hormone metabolic process",
                            "branched-chain amino acid metabolic process"))+
  scale_fill_manual(values=c("#EC0000","#0055BE"),name = "",labels=c("Upregulated","Downregulated"))+
  geom_line(size=4,color="lightgrey")+
  geom_pie_glyph(slices = c('deg_up_num', 'deg_down_num'),radius=0.4)+
  geom_point(size=5,colour="white")+
  theme_bw()+
  xlab("")+
  ylab("")+
  labs(tag="B")+
  ggtitle("GO Biological Process")+
  scale_x_discrete(limits=c("rearing","testing"),labels=c("Developmental 
Temp.","Testing 
Temp."))+
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title= element_text(size=10),
        plot.title = element_text(size = 11, face = "bold"),
        #        legend.margin=margin(t = 0, unit='cm'),
        legend.position = "right")
#        legend.spacing.y = unit(0.5, 'cm'))
bpplot_legend

legend <- get_legend(bpplot_legend)



### now the same for KEGG results 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_Enrichment_analysis_outputfiles")
kegg_testing <- read.csv("testingtempKEGG.csv")
kegg_rearing <- read.csv("rearingtempKEGG.csv")

kegg_testing$treatment <- "testing"
kegg_rearing$treatment <- "rearing"
dfboth_kegg <- rbind(kegg_rearing,kegg_testing)

sort.dfboth_kegg <- with(dfboth_kegg,  dfboth_kegg[order(Description) , ])

linegroup_kegg <- c("UDP2","UDP1","circ2","UDP2","UDP2","UDP1","HSP2","HSP1",
                    "UDP2","UDP1","UDP2","UDP1","UDP2","UDP1","HSP2","HSP1",
                    "UDP2","UDP1")
sort.dfboth_kegg$linegroup_kegg <- linegroup_kegg

keggplot<-ggplot(data=sort.dfboth_kegg,aes(x=treatment,y=Description,group=linegroup_kegg))+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25),
                   limits=c("Protein processing in endoplasmic reticulum - Pieris rapae (cabbage white)",
                            "Longevity regulating pathway - multiple species - Pieris rapae (cabbage white)",
                            "Retinol metabolism - Pieris rapae (cabbage white)",
                            "Pentose and glucuronate interconversions - Pieris rapae (cabbage white)",
                            "Metabolism of xenobiotics by cytochrome P450 - Pieris rapae (cabbage white)",
                            "Drug metabolism - other enzymes - Pieris rapae (cabbage white)",
                            "Ascorbate and aldarate metabolism - Pieris rapae (cabbage white)",
                            "Porphyrin metabolism - Pieris rapae (cabbage white)",
                            "Circadian rhythm - fly - Pieris rapae (cabbage white)"))+
  scale_fill_manual(values=c("#EC0000","#0055BE"),name = "",labels=c("Upregulated","Downregulated"))+
  geom_line(size=4,color="lightgrey")+
  geom_pie_glyph(slices = c('deg_up_num', 'deg_down_num'),radius=0.4)+
  geom_point(size=5,colour="white")+
  theme_bw()+
  xlab("")+
  ylab("")+
 # labs(tag="C")+
  ggtitle("(C) KEGG Pathway")+
  scale_x_discrete(limits=c("rearing","testing"),labels=c("Dev. 
Temp.","Short-term Acc. 
Temp."))+
  geom_label(aes(label=Count),fill=NA,label.size=NA,size=3)+
  #  scale_size_continuous(breaks = c(3, 6, 9, 12))
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title= element_text(size=10),
        plot.title = element_text(size = 11, face = "bold"),
        #        legend.margin=margin(t = 0, unit='cm'),
        legend.position = "")
#        legend.spacing.y = unit(0.5, 'cm'))
keggplot

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("enrichment_results_donut.png",plot=grid.arrange(mfplot,bpplot,keggplot, legend, ncol=4,nrow=1, widths=c(0.3,0.3,0.3,0.1),heights=c(1)),dpi=600,
       units='in',width=13,height=8)



