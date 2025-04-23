


## use clustering to look at the interaction genes
# https://tavareshugo.github.io/data-carpentry-rnaseq/04a_explore_test_results.html
# https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

## get list of candidate genes
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
interaction_sig <- read.csv("2factordesign_interaction_sig.gene.csv")
norm_counts <- read.csv("counts_transformed.csv")
sample_info <- read.csv("sample_info.csv")

## summarize counts
norm_counts_mean <- norm_counts %>% 
  # convert to long format
  pivot_longer(cols = starts_with("sample_"), names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # filter to retain only genes of interest
  filter(gene %in% interaction_sig$gene) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(gene, condition) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()

head(norm_counts_mean)

# Create a matrix
hclust_matrix <- norm_counts %>% 
  select(-gene) %>% 
  as.matrix()
hclust_matrix
t(hclust_matrix)
# assign rownames
rownames(hclust_matrix) <- norm_counts$gene
hclust_matrix <- hclust_matrix[interaction_sig$gene, ]
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scaling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()
hclust_matrix
gene_dist <- dist(hclust_matrix)
gene_dist
gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 5, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

## https://www.stat.cmu.edu/~ryantibs/datamining/lectures/06-clus3.pdf => advice on how to pick the number of clusters
install.packages("cluster")
library("cluster")
# https://r-tastic.co.uk/post/optimal-number-of-clusters/ 
install.packages("knitr","gclus","dpylr")
library(dplyr)
library(knitr)
library(gclus)

### elbow method of picking cluster number => As a rule of thumb, you pick the 
# number for which you see a significant decrease in the within-cluster dissimilarity,
# or so called ‘elbow’
wss <- (nrow(hclust_matrix)-1)*sum(apply(hclust_matrix,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(hclust_matrix,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
# I see moderate decreases up to about 8 clusters

## version with 5 clusters
gene_cluster_k5 <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
gene_cluster_k5 <- merge(gene_cluster_k5,interaction_sig,by="gene")
gene_cluster_k5 <- merge(gene_cluster_k5,interaction_sig,by="gene")

## version with 8 clusters
gene_cluster_k8 <- cutree(gene_hclust, k = 8) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
gene_cluster_k8 <- merge(gene_cluster_k8,interaction_sig,by="gene")
View(gene_cluster_k8)
gene_cluster_k8 %>% count(cluster) ## tells the number of genes in each cluster

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/RNA-Seq workflow/R output")
write.csv(gene_cluster_k5,"interactiongene_clusters_k5.csv",row.names = F)
write.csv(gene_cluster_k8,"interactiongene_clusters_k8.csv",row.names = F)

norm_counts_cluster_k5 <- norm_counts_mean %>% 
  inner_join(gene_cluster_k5, by = "gene")

norm_counts_cluster_k8 <- norm_counts_mean %>% 
  inner_join(gene_cluster_k8, by = "gene")


## adds rearing temp and testing temp columns 
Devtemp <- substr(norm_counts_cluster_k8$condition, 1, 2) # this takes the first 2 characters
temps <- c("16","30","30","16") 
N <- length(norm_counts_cluster_k8$gene)/4
Testingtemp <- rep(temps, times=N)
norm_counts_cluster_k8$Testingtemp <- Testingtemp
norm_counts_cluster_k8$Devtemp <- Devtemp
View(norm_counts_cluster_k8)

### re-naming the facets titles
cluster.labsk5 <- c("Cluster 1","Cluster 2","Cluster 3",
                    "Cluster 4","Cluster 5")
names(cluster.labsk5) <- c(1,2,3,4,5)

cluster.labsk8 <- c("Cluster 1","Cluster 2","Cluster 3",
                    "Cluster 4","Cluster 5", "Cluster 6","Cluster 7",
                    "Cluster 8")
names(cluster.labsk8) <- c(1,2,3,4,5,6,7,8)

View(norm_counts_cluster_k8)
norms_counts_cluster_k8_subset <- norm_counts_cluster_k8[,c(1,6,14,13,3,5)]
View(norms_counts_cluster_k8_subset)

## save this normalized mean scale count file 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Script_outputfiles/Prapae_IntGene_Clustering_outputfiles")
write.csv(norms_counts_cluster_k8_subset,"interactiongene_clusters_k8_MeanScaledCounts.csv",row.names = F)


## new clusterplot in reaction norm -style with faint lines for each gene and thick lines for averaging across genes each cluster
clusterplot_k5 <- norm_counts_cluster_k5 %>% 
  ggplot(aes(Testingtemp, mean_cts_scaled, color = Devtemp,group=Devtemp)) +
  stat_summary(fun.y=mean,size=0.5,aes(group=Devtemp))+
  #geom_point(size=2,alpha=0.3,position=position_jitter(w = 0.2,h = 0)) +
  stat_summary(fun.y=mean,geom='line',aes(group=factor(Devtemp),colour=Devtemp),size=0.7)+
  stat_summary(fun.y=mean,geom='line',aes(group=interaction(gene,Devtemp),colour=Devtemp),size=0.35,alpha=0.25)+
  stat_summary(fun.data = mean_se,geom="errorbar", width = 0.2,size=0.7) +
  facet_grid(cols = vars(cluster),labeller = labeller(cluster = cluster.labsk5))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title= element_text(size=11),
        plot.title = element_text(size = 11, face = "bold"))+
  theme(strip.text.x = element_text(size = 11,face="bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0, linetype="solid"))+
  scale_color_manual(values = c("darkblue", "darkred"), 
                     name = "Rearing 
Temperature (°C)")+
  xlab("Testing Temperature (°C)")+
  ylab("Mean Counts Scaled")

clusterplot_k5
View(norm_counts_cluster_k8)
clusterplot_k8 <- norm_counts_cluster_k8 %>% 
  ggplot(aes(Testingtemp, mean_cts_scaled, color = Devtemp,group=Devtemp)) +
  stat_summary(fun.data = mean_se,geom="errorbar", fun.args = list(mult = 2), width = 0.2,size=0.7) +
  stat_summary(fun.y=mean,size=0.5,aes(group=Devtemp,shape=Devtemp))+
  stat_summary(fun.y=mean,geom='line',aes(group=factor(Devtemp),colour=Devtemp,linetype=Devtemp),size=0.7)+
  stat_summary(fun.y=mean,geom='line',aes(group=interaction(gene,Devtemp),colour=Devtemp,linetype=Devtemp),size=0.35,alpha=0.25)+
  facet_grid(cols = vars(cluster),labeller = labeller(cluster = cluster.labsk8))+
  theme_bw()+
  scale_x_discrete(breaks = c("16","30"),
                   labels = c("16±5","30±5"))+
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title= element_text(size=11),
        plot.title = element_text(size = 11, face = "bold"))+
  theme(strip.text.x = element_text(size = 11,face="bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0, linetype="solid"))+
  scale_color_manual(values = c("darkblue", "darkred"),
                    breaks = c("16","30"),
                     labels = c("16±5","30±5"),
                     name = "Developmental 
Temperature (°C)")+
  scale_shape_manual(values = c(19, 21),
                     breaks = c("16","30"),
                     labels = c("16±5","30±5"),
                     name = "Developmental 
Temperature (°C)")+
  scale_linetype_manual(values = c("solid", "dashed"),
                        breaks = c("16","30"),
                        labels = c("16±5","30±5"),
                        name = "Developmental 
Temperature (°C)")+
  xlab("Testing Temperature (°C)")+
  ylab("Mean Counts Scaled")+
  guides(color = guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 1,direction="horizontal",override.aes = list(size=0.75)),
         shape = guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 1,direction="horizontal"),
         linetype=guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 3,direction="horizontal",override.aes = list(size=0.75)))

clusterplot_k8

clusterplot_k8_cluster3 <- norm_counts_cluster_k8[norm_counts_cluster_k8$cluster==3,] %>% 
  ggplot(aes(Testingtemp, mean_cts_scaled, color = Devtemp,group=Devtemp)) +
  stat_summary(fun.data = mean_se,geom="errorbar", fun.args = list(mult = 2), width = 0.1,size=0.5) +
  stat_summary(fun.y=mean,size=0.4,aes(group=Devtemp,shape=Devtemp))+
  stat_summary(fun.y=mean,geom='line',aes(group=factor(Devtemp),colour=Devtemp,linetype=Devtemp),size=0.5)+
  stat_summary(fun.y=mean,geom='line',aes(group=interaction(gene,Devtemp),colour=Devtemp,linetype=Devtemp),size=0.35,alpha=0.25)+
  facet_grid(cols = vars(cluster),labeller = labeller(cluster = cluster.labsk8))+
  theme_bw()+
  scale_x_discrete(breaks = c("16","30"),
                   labels = c("16±5","30±5"))+
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title= element_text(size=9),
        plot.title = element_text(size = 10, face = "bold"))+
  theme(strip.text.x = element_text(size = 10,face="bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0, linetype="solid"))+
  scale_color_manual(values = c("darkblue", "darkred"),
                     breaks = c("16","30"),
                     labels = c("16±5","30±5"),
                     name = "Developmental 
Temperature (°C)")+
  scale_shape_manual(values = c(19, 21),
                     breaks = c("16","30"),
                     labels = c("16±5","30±5"),
                     name = "Developmental 
Temperature (°C)")+
  scale_linetype_manual(values = c("solid", "dashed"),
                        breaks = c("16","30"),
                        labels = c("16±5","30±5"),
                        name = "Developmental 
Temperature (°C)")+
  xlab("Testing Temperature (°C)")+
  ylab("Mean Counts Scaled")+
  guides(color = guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 1,direction="horizontal",override.aes = list(size=0.5)),
         shape = guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 1,direction="horizontal"),
         linetype=guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 3,direction="horizontal",override.aes = list(size=0.5)))+
  theme(legend.position = "top")

clusterplot_k8_cluster3


clusterplot_k8_verticalclusters <- norm_counts_cluster_k8 %>% 
  ggplot(aes(Testingtemp, mean_cts_scaled, color = Devtemp,group=Devtemp)) +
  stat_summary(fun.data = mean_se,geom="errorbar", fun.args = list(mult = 2), width = 0.1,size=0.5) +
  stat_summary(fun.y=mean,size=0.3,aes(group=Devtemp,shape=Devtemp))+
  stat_summary(fun.y=mean,geom='line',aes(group=factor(Devtemp),colour=Devtemp,linetype=Devtemp),size=0.5)+
  stat_summary(fun.y=mean,geom='line',aes(group=interaction(gene,Devtemp),colour=Devtemp,linetype=Devtemp),size=0.25,alpha=0.25)+
  facet_grid(rows = vars(cluster),labeller = labeller(cluster = cluster.labsk8))+
  theme_bw()+
  scale_x_discrete(breaks = c("16","30"),
                   labels = c("16±5","30±5"))+
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size=9),
        legend.text = element_text(size=8),
        legend.title= element_text(size=9),
        plot.title = element_text(size = 9, face = "bold"))+
  theme(strip.text.x = element_text(size = 10,face="bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0, linetype="solid"))+
  scale_color_manual(values = c("darkblue", "darkred"),
                     breaks = c("16","30"),
                     labels = c("16±5","30±5"),
                     name = "Developmental 
Temperature (°C)")+
  scale_shape_manual(values = c(19, 21),
                     breaks = c("16","30"),
                     labels = c("16±5","30±5"),
                     name = "Developmental 
Temperature (°C)")+
  scale_linetype_manual(values = c("solid", "dashed"),
                        breaks = c("16","30"),
                        labels = c("16±5","30±5"),
                        name = "Developmental 
Temperature (°C)")+
  xlab("Testing Temperature (°C)")+
  ylab("Mean Counts Scaled")+
  guides(color = guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 1,direction="horizontal",override.aes = list(size=0.75)),
         shape = guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 1,direction="horizontal"),
         linetype=guide_legend(title = "Developmental 
Temperature (°C)", order = 1, ncol = 3,direction="horizontal",override.aes = list(size=0.75)))+
  theme(legend.position = "top")

clusterplot_k8_verticalclusters


setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/INTBio grant work/Writing/Final Scripts/Final Figures")
ggsave("clusterplot_k8_reactionnorm.png",plot=clusterplot_k8,dpi=600,
       units='in',width=9,height=3)

ggsave("clusterplot_k8_reactionnorm_vertical.png",plot=clusterplot_k8_verticalclusters,dpi=600,
       units='in',width=2.5,height=11)

ggsave("clusterplot_k8_reactionnorm_cluster3.png",plot=clusterplot_k8_cluster3,dpi=600,
       units='in',width=4,height=4)



