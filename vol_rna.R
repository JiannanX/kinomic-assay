
library(ggplot2)
library(clusterProfiler)
library(edgeR)
library(dplyr)
library(ggrepel)

# setwd("/Users/nicole/Desktop/11.6/phe&asp")
setwd("D:/WeChat/WeChat Files/wxid_p8frykgnz4m522/FileStorage/File/2023-01")

ptk = read.csv("uniprot_name_ptk.csv", header = T)[,c(2,3)]
stk = read.csv("uniprot_name_stk.csv", header = T)[,c(2,3)]

colnames(stk) = colnames(ptk) = c("kinase", "kinase_ID")
all_kinase = rbind(ptk,stk)

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
gene_all <- bitr(all_kinase[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)
gene2_all = merge(gene_all,all_kinase, by.x = 'UNIPROT', by.y = "kinase_ID")


# all the kinase significant in ptk.stk, ptk vs stk in ptk&ast- all

# setwd("/Users/nicole/Desktop/KS/kinase function")
kinases_name <- read.csv2("ptk&stk-all.csv", header = F)

match_result = bitr(gene2_all$ENTREZID,fromType = 'ENTREZID', toType = 'ENSEMBL',OrgDb = GO_database)
match_result = merge(match_result, gene2_all, by.x = "ENTREZID", by.y = "ENTREZID")
match_result = match_result[,-3]

match_result2 = merge(kinases_name, gene2_all, by.x = colnames(kinases_name)[1], by.y = "kinase")
match_result2 <- bitr(match_result2$ENTREZID,fromType = 'ENTREZID', toType = 'ENSEMBL',OrgDb = GO_database)
match_result2 = merge(match_result2, gene2_all, by.x = "ENTREZID", by.y = "ENTREZID")
match_result2 = match_result2[,-3]
colnames(match_result) = colnames(match_result2) = c("ENTREZID","ENSEMBL","V3")

# rna-seq for T2D VS non-diabetic

counts <- read.csv("T2D&ND.csv", header = T)
rownames(counts)=counts[,1]
counts=counts[,-1]
group <- c(rep(1,18),rep(2,39))
y = DGEList(counts=counts, group=group)

keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
et <- topTags(et, n=100000)
et <- as.data.frame(et)
et <- cbind(rownames(et),et)

colnames(et) <- c("gene_id", "log2FoldChange", "log2CPM", "PValue", "FDR")

write.csv(et, "DEG-all.gene.csv")

etSig <- et[which(et$PValue < 0.05 & abs(et$log2FoldChange) > 1),]

etSig[which(etSig$log2FoldChange > 0), "up_down"] <- "Up"
etSig[which(etSig$log2FoldChange < 0), "up_down"] <- "Down"

write.csv(etSig, "Pval0.05-FC2-edgeR.csv")

# volcano plot

v1<-read.csv("DEG-all.gene.csv", header = T)
library(stringr)
v1$gene_id <- str_replace(v1$gene_id, pattern = ".[0-9]+$",replacement = "")
v1 = v1[,-1]

v3 <- v1 %>% right_join(match_result2,by=c("gene_id"="ENSEMBL")) %>% na.omit
v3$change<-as.factor(ifelse(v3$PValue<0.05 & abs(v3$log2FoldChange)>1, ifelse(v3$log2FoldChange>1,"up","down"),"not significant"))


v2 = v1 %>% right_join(match_result,by=c("gene_id"="ENSEMBL")) %>% na.omit
v2$change<-as.factor(ifelse(v2$PValue<0.05 , "significant","not significant"))
v2[!(v2$V3 %in% v3$V3),7] = NA
v2[v2$change == "not significant",7] = NA

# no fold change filter

options(ggrepel.max.overlaps = Inf)
tiff("volcano_FC1_1.tiff",width = 1500, height = 1800, res = 200)
ggplot(data = v2, aes(x = log2FoldChange, y = -log10(PValue), color = change)) +
  geom_point(alpha=0.7, size = 2) + xlim(-2,2) + ylim(0,4) + 
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_hline(yintercept=-log10(0.05) ,linetype=4) +
  scale_color_manual(name = "", values = c("red", "grey"), 
                     limits = c("significant", "not significant")) +
  geom_label_repel(aes(label = v2$V3), 
                   fontface="bold", color="black", min.segment.length = 0, 
                   segment.colour = "black")
dev.off()

# Foldchange > 2 filter

v2$change<-as.factor(ifelse(v2$PValue<0.05 & abs(v2$log2FoldChange)>1, ifelse(v2$log2FoldChange>1,"up","down"),"not significant"))
v3[v3$change == "not significant",7] = NA
v2[!(v2$V3 %in% v3$V3),7] = NA

options(ggrepel.max.overlaps = Inf)
tiff("volcano_FC1_2.tiff",width = 1500, height = 1800, res = 200)
ggplot(data = v2, aes(x = log2FoldChange, y = -log10(PValue), color = change)) +
  geom_point(alpha=0.7, size = 2) + xlim(-2,2) + ylim(0,4) + 
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_hline(yintercept=-log10(0.05) ,linetype=4) +
  geom_vline(xintercept=c(-1,1) ,linetype=4) +
  scale_color_manual(name = "", values = c("red", "blue", "grey"), 
                     limits = c("up", "down", "not significant")) +
  geom_label_repel(aes(label = v2$V3), 
      fontface="bold", color="black", min.segment.length = 0, 
 segment.colour = "black")
dev.off()
