library(clusterProfiler)
library(ggplot2)
library(ReactomePA)

## Integration of kinases that is significantly changed in either Ex4, Phe or Asp, ptk 10mins / 120mins and stk 10mins / 120mins.

setwd("/Users/nicole/Desktop/11.6/ptk")
ptk10 = read.csv("ptk10all.csv", header = T)
ptk120 = read.csv("ptk120all.csv", header = T)

ptk = union(as.character(ptk10[,2]), as.character(ptk120[,2]))
uniprot_ptk = read.csv("uniprot_name_ptk.csv", header = T)
ptk = uniprot_ptk[uniprot_ptk[,2] %in% ptk,c(2,3)]

setwd("/Users/nicole/Desktop/11.6/stk")
stk10 = read.csv("stk10all.csv", header = T)
stk120 = read.csv("stk120all.csv", header = T)

stk = union(as.character(stk10[,2]), as.character(stk120[,2]))
uniprot_stk = read.csv("uniprot_name_stk.csv", header = T)
stk = uniprot_stk[uniprot_stk[,2] %in% stk,c(2,3)]

colnames(stk) = c("kinase", "kinase_ID")
all_kinase = rbind(ptk,stk)


## Retrieve kinase ENTREZID

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

gene_ptk <- bitr(ptk[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)
gene_stk <- bitr(stk[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)
gene_all <- bitr(all_kinase[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)

gene2_ptk = merge(gene_ptk,ptk, by.x = 'UNIPROT', by.y = "kinase_ID")
gene2_stk = merge(gene_stk,stk, by.x = 'UNIPROT', by.y = "kinase_ID")
gene2_all = merge(gene_all,all_kinase, by.x = 'UNIPROT', by.y = "kinase_ID")


## GO barplot for ptk / stk and their both - Biological Process / Cellular Components 

BP_ptk <- enrichGO(gene2_ptk$ENTREZID,OrgDb = GO_database,ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(BP_ptk)
tiff("BP_ptk.tiff",width = 1500, height = 1800, res = 200)
mutate(BP_ptk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

BP_stk <- enrichGO(gene2_stk$ENTREZID,OrgDb = GO_database,ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(BP_stk)
tiff("BP_stk.tiff",width = 1500, height = 1800, res = 200)
mutate(BP_stk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

BP_all <- enrichGO(gene2_all$ENTREZID,OrgDb = GO_database,ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(BP_all)
tiff("BP_all.tiff",width = 1500, height = 1800, res = 200)
mutate(BP_all, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()


CC_ptk <- enrichGO(gene2_ptk$ENTREZID,OrgDb = GO_database,ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(CC_ptk)
tiff("CC_ptk.tiff",width = 1500, height = 1800, res = 200)
mutate(CC_ptk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

CC_stk <- enrichGO(gene2_stk$ENTREZID,OrgDb = GO_database,ont = "CC", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(CC_stk)
tiff("CC_stk.tiff",width = 1500, height = 1800, res = 200)
mutate(CC_stk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

CC_all <- enrichGO(gene2_all$ENTREZID,OrgDb = GO_database,ont = "CC", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(CC_all)
tiff("CC_all.tiff",width = 1500, height = 1800, res = 200)
mutate(CC_all, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

## KEGG dotplot for ptk / stk and their both

KEGG_ptk<-enrichKEGG(gene2_ptk$ENTREZID, organism = KEGG_database , pvalueCutoff = 0.05)
head(KEGG_ptk,10)
tiff("KEGG_ptk_ratio.tiff",width = 1500, height = 1800, res = 200)
dotplot(KEGG_ptk, showCategory=10)
dev.off()

tiff("KEGG_ptk_qscore.tiff",width = 1500, height = 1800, res = 200)
mutate(KEGG_ptk, qscore = -log(p.adjust, base=10)) %>% 
  dotplot(showCategory=10, x="qscore",font.size = 12)
dev.off()


KEGG_stk<-enrichKEGG(gene2_stk$ENTREZID, organism = KEGG_database, pvalueCutoff = 0.05)
head(KEGG_stk,10)
tiff("KEGG_stk_ratio.tiff",width = 1500, height = 1800, res = 200)
dotplot(KEGG_stk, showCategory=10)
dev.off()

tiff("KEGG_stk_qscore.tiff",width = 1500, height = 1800, res = 200)
mutate(KEGG_stk, qscore = -log(p.adjust, base=10)) %>% 
  dotplot(showCategory=10, x="qscore")
dev.off()


KEGG_all<-enrichKEGG(gene2_all$ENTREZID, organism = KEGG_database, pvalueCutoff = 0.05)
head(KEGG_all)
tiff("KEGG_all_ratio.tiff",width = 1500, height = 1800, res = 200)
dotplot(KEGG_all, showCategory=10)
dev.off()

tiff("KEGG_all_qscore.tiff",width = 1500, height = 1800, res = 200)
mutate(KEGG_all, qscore = -log(p.adjust, base=10)) %>% 
  dotplot(showCategory=10, x="qscore")
dev.off()

## Reactome barplots for ptk / stk and their both


RE_ptk<-enrichPathway(gene2_ptk$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(RE_ptk)
tiff("RE_ptk.tiff",width = 1500, height = 1800, res = 200)
mutate(RE_ptk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

RE_stk<-enrichPathway(gene2_stk$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(RE_stk)
tiff("RE_stk.tiff",width = 1500, height = 1800, res = 200)
mutate(RE_stk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

RE_all<-enrichPathway(gene2_all$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(RE_all)
tiff("RE_all.tiff",width = 1500, height = 1800, res = 200)
mutate(RE_all, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()


## VENN plot 1

head(KEGG_ptk,10)

ptk_select_1 = strsplit("3643/1436/3480/5159/2260/2321/2263/2264/2261/3791/2322/5599/5601/5602/5604/4915", "/")[1]
ptk_select_2 = strsplit("9252/3643/1436/3480/5159/2260/2321/2263/2264/2261/3791/2322/5599/5601/5602/5604/4915", "/")[1]
ptk_select_3 = strsplit("3717/3643/1436/3480/5159/2260/2321/2263/2264/2261/3791/2322/3718/5604/4915", "/")[1]
ptk_select_4 = strsplit("3717/3480/5159/2263/2261/3791/5604", "/")[1]
ptk_select_5 = strsplit("5979/5159/2260/2321/2263/2264/2261/3791/4915", "/")[1]

ptk_select_1 = as.data.frame(ptk_select_1)
ptk_select_2 = as.data.frame(ptk_select_2)
ptk_select_3 = as.data.frame(ptk_select_3)
ptk_select_4 = as.data.frame(ptk_select_4)
ptk_select_5 = as.data.frame(ptk_select_5)

colnames(ptk_select_1) = colnames(ptk_select_2) = colnames(ptk_select_3) = colnames(ptk_select_4) = colnames(ptk_select_5) = "select"
ptk_select_1 = merge(ptk_select_1,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_2 = merge(ptk_select_2,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_3 = merge(ptk_select_3,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_4 = merge(ptk_select_4,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_5 = merge(ptk_select_5,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]


head(KEGG_stk,10)

stk_select_1 = strsplit("3551/5603/5170/8767/9252/5894/673/5595/5594/207/208/5599/5601/9261/2932/6197/6300/5602/5598/816/818/817/6196/6195/5600/1432/27330/815/10000", "/")[1]
stk_select_2 = strsplit("6446/3551/1147/5603/5170/5894/369/673/1017/5595/5594/207/208/5599/5601/1454/102800317/5347/6300/5602/5562/5600/6794/1432/23678/100533105/10110/10000", "/")[1]
stk_select_3 = strsplit("10298/5894/5579/369/673/5578/6198/5595/5594/207/208/2475/5599/5601/2932/5602/816/818/817/57144/6199/815/10000", "/")[1]
stk_select_4 = strsplit("3551/1147/5603/9252/8986/5894/5579/369/673/5578/5566/5568/5567/5595/5594/207/208/5599/5601/5606/9261/6197/5608/6300/5602/5598/6196/6195/5600/1432/7867/8491/8569/27330/10000", "/")[1]
stk_select_5 = strsplit("5603/5894/369/1019/673/5566/5568/5567/6198/5595/5594/207/208/2475/5599/5601/6300/5602/5600/1432/6199/10000", "/")[1]

stk_select_1 = as.data.frame(stk_select_1)
stk_select_2 = as.data.frame(stk_select_2)
stk_select_3 = as.data.frame(stk_select_3)
stk_select_4 = as.data.frame(stk_select_4)
stk_select_5 = as.data.frame(stk_select_5)

colnames(stk_select_1) = colnames(stk_select_2) = colnames(stk_select_3) = colnames(stk_select_4) = colnames(stk_select_5) = "select"

stk_select_1 = merge(stk_select_1,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_2 = merge(stk_select_2,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_3 = merge(stk_select_3,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_4 = merge(stk_select_4,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_5 = merge(stk_select_5,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]


# VENN plot2

ptk_select_1 = strsplit("3643/1436/3480/5159/2260/2321/2263/2264/2261/3791/2322/5599/5601/5602/5604/4915", "/")[1]
ptk_select_2 = strsplit("9252/3643/1436/3480/5159/2260/2321/2263/2264/2261/3791/2322/5599/5601/5602/5604/4915", "/")[1]
ptk_select_3 = strsplit("3717/3643/1436/3480/5159/2260/2321/2263/2264/2261/3791/2322/3718/5604/4915", "/")[1]

ptk_select_1 = as.data.frame(ptk_select_1)
ptk_select_2 = as.data.frame(ptk_select_2)
ptk_select_3 = as.data.frame(ptk_select_3)

colnames(ptk_select_1) = colnames(ptk_select_2) = colnames(ptk_select_3) = "select"

ptk_select_1 = merge(ptk_select_1,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_2 = merge(ptk_select_2,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_3 = merge(ptk_select_3,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]


stk_select_1 = strsplit("3551/1147/5603/9252/8986/5894/5579/369/673/5578/5566/5568/5567/5595/5594/207/208/5599/5601/5606/9261/6197/5608/6300/5602/5598/6196/6195/5600/1432/7867/8491/8569/27330/10000", "/")[1]
stk_select_2 = strsplit("3551/5603/5170/8767/9252/5894/673/5595/5594/207/208/5599/5601/9261/2932/6197/6300/5602/5598/816/818/817/6196/6195/5600/1432/27330/815/10000", "/")[1]
stk_select_3 = strsplit("6446/3551/1147/5603/5170/5894/369/673/1017/5595/5594/207/208/5599/5601/1454/102800317/5347/6300/5602/5562/5600/6794/1432/23678/100533105/10110/10000", "/")[1]

stk_select_1 = as.data.frame(stk_select_1)
stk_select_2 = as.data.frame(stk_select_2)
stk_select_3 = as.data.frame(stk_select_3)

colnames(stk_select_1) = colnames(stk_select_2) = colnames(stk_select_3) = "select"

stk_select_1 = merge(stk_select_1,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_2 = merge(stk_select_2,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_3 = merge(stk_select_3,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]


## Integration of kinases that is significantly changed in either Phe or Asp normalized to Ex4, ptk 10mins / 120mins and stk 10mins / 120mins.

setwd("/Users/nicole/Desktop/11.6/phe&asp")
ptk10 = read.csv("ptk10all.csv", header = T)
ptk120 = read.csv("ptk120all.csv", header = T)

ptk = union(as.character(ptk10[,2]), as.character(ptk120[,2]))
uniprot_ptk = read.csv("uniprot_name_ptk.csv", header = T)
ptk = uniprot_ptk[uniprot_ptk[,2] %in% ptk,c(2,3)]

stk10 = read.csv("stk10all.csv", header = T)
stk120 = read.csv("stk120all.csv", header = T)

stk = union(as.character(stk10[,2]), as.character(stk120[,2]))
uniprot_stk = read.csv("uniprot_name_stk.csv", header = T)
stk = uniprot_stk[uniprot_stk[,2] %in% stk,c(2,3)]

colnames(stk) = c("kinase", "kinase_ID")
all_kinase = rbind(ptk,stk)


GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'

gene_ptk <- bitr(ptk[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)
gene_stk <- bitr(stk[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)
gene_all <- bitr(all_kinase[,2],fromType = 'UNIPROT',toType = 'ENTREZID',OrgDb = GO_database)

gene2_ptk = merge(gene_ptk,ptk, by.x = 'UNIPROT', by.y = "kinase_ID")
gene2_stk = merge(gene_stk,stk, by.x = 'UNIPROT', by.y = "kinase_ID")
gene2_all = merge(gene_all,all_kinase, by.x = 'UNIPROT', by.y = "kinase_ID")


## GO barplot for ptk / stk and their both - Biological Process / Cellular Components 

BP_ptk <- enrichGO(gene2_ptk$ENTREZID,OrgDb = GO_database,ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(BP_ptk)
tiff("BP_ptk.tiff",width = 1500, height = 1800, res = 200)
mutate(BP_ptk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

BP_stk <- enrichGO(gene2_stk$ENTREZID,OrgDb = GO_database,ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(BP_stk)
tiff("BP_stk.tiff",width = 1500, height = 1800, res = 200)
mutate(BP_stk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

BP_all <- enrichGO(gene2_all$ENTREZID,OrgDb = GO_database,ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(BP_all)
tiff("BP_all.tiff",width = 1500, height = 1800, res = 200)
mutate(BP_all, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()


CC_ptk <- enrichGO(gene2_ptk$ENTREZID,OrgDb = GO_database,ont = "CC", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(CC_ptk)
tiff("CC_ptk.tiff",width = 1500, height = 1800, res = 200)
mutate(CC_ptk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

CC_stk <- enrichGO(gene2_stk$ENTREZID,OrgDb = GO_database,ont = "CC", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(CC_stk)
tiff("CC_stk.tiff",width = 1500, height = 1800, res = 200)
mutate(CC_stk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

CC_all <- enrichGO(gene2_all$ENTREZID,OrgDb = GO_database,ont = "CC", pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable = TRUE)
head(CC_all)
tiff("CC_all.tiff",width = 1500, height = 1800, res = 200)
mutate(CC_all, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()


## KEGG dotplot for ptk / stk and their both

KEGG_ptk<-enrichKEGG(gene2_ptk$ENTREZID, organism = KEGG_database , pvalueCutoff = 0.05)
head(KEGG_ptk)
tiff("KEGG_ptk_ratio.tiff",width = 1500, height = 1800, res = 200)
dotplot(KEGG_ptk, showCategory=10)
dev.off()

tiff("KEGG_ptk_qscore.tiff",width = 1500, height = 1800, res = 200)
mutate(KEGG_ptk, qscore = -log(p.adjust, base=10)) %>% 
  dotplot(showCategory=10, x="qscore")
dev.off()

KEGG_stk<-enrichKEGG(gene2_stk$ENTREZID, organism = KEGG_database, pvalueCutoff = 0.05)
head(KEGG_stk,10)
tiff("KEGG_stk_ratio.tiff",width = 1500, height = 1800, res = 200)
dotplot(KEGG_stk, showCategory=10)
dev.off()

tiff("KEGG_stk_qscore.tiff",width = 1500, height = 1800, res = 200)
mutate(KEGG_stk, qscore = -log(p.adjust, base=10)) %>% 
  dotplot(showCategory=10, x="qscore")
dev.off()

KEGG_all<-enrichKEGG(gene2_all$ENTREZID, organism = KEGG_database, pvalueCutoff = 0.05)
head(KEGG_all)
tiff("KEGG_all_ratio.tiff",width = 1500, height = 1800, res = 200)
dotplot(KEGG_all, showCategory=10)
dev.off()

tiff("KEGG_all_qscore.tiff",width = 1500, height = 1800, res = 200)
mutate(KEGG_all, qscore = -log(p.adjust, base=10)) %>% 
  dotplot(showCategory=10, x="qscore")
dev.off()


## Reactome barplot for ptk / stk and their both

RE_ptk<-enrichPathway(gene2_ptk$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(RE_ptk)
tiff("RE_ptk.tiff",width = 1500, height = 1800, res = 200)
mutate(RE_ptk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

RE_stk<-enrichPathway(gene2_stk$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(RE_stk)
tiff("RE_stk.tiff",width = 1500, height = 1800, res = 200)
mutate(RE_stk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()

RE_all<-enrichPathway(gene2_all$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(RE_all)
tiff("RE_all.tiff",width = 1500, height = 1800, res = 200)
mutate(RE_all, qscore = -log(p.adjust, base=10)) %>% 
  barplot(showCategory=10, x="qscore")
dev.off()


# VENN 1

KEGG_ptk<-enrichKEGG(gene2_ptk$ENTREZID, organism = KEGG_database , pvalueCutoff = 0.05)
head(KEGG_ptk,10)

ptk_select_1 = strsplit("3643/1436/3480/5159/2260/5156/2321/2263/2261/1969/3791/2322/5599/5601/5602", "/")[1]
ptk_select_2 = strsplit("3643/1436/3480/5159/2260/5156/2321/2263/2261/1969/3791/2322/5599/5601/5602", "/")[1]
ptk_select_3 = strsplit("3717/3643/1436/3480/5159/2260/5156/2321/2263/2261/1969/3791/2322/3718", "/")[1]

ptk_select_1 = as.data.frame(ptk_select_1)
ptk_select_2 = as.data.frame(ptk_select_2)
ptk_select_3 = as.data.frame(ptk_select_3)

colnames(ptk_select_1) = colnames(ptk_select_2) = colnames(ptk_select_3) = "select"
ptk_select_1 = merge(ptk_select_1,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_2 = merge(ptk_select_2,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_3 = merge(ptk_select_3,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]


KEGG_stk<-enrichKEGG(gene2_stk$ENTREZID, organism = KEGG_database , pvalueCutoff = 0.05)
head(KEGG_stk,20)

stk_select_1 = strsplit("10298/5894/5582/5579/369/673/5578/2932/816/57144", "/")[1]
stk_select_2 = strsplit("3551/1147/9252/5894/5582/5579/369/673/5578/5606/5608/7867/8491/51701", "/")[1]

stk_select_1 = as.data.frame(stk_select_1)
stk_select_2 = as.data.frame(stk_select_2)

colnames(stk_select_1) = colnames(stk_select_2) = "select"

stk_select_1 = merge(stk_select_1,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
stk_select_2 = merge(stk_select_2,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]


# VENN 2

ptk_select_1 = strsplit("3643/1436/3480/5159/2260/5156/2321/2263/2261/1969/3791/2322/5599/5601/5602", "/")[1]
ptk_select_2 = strsplit("3643/1436/3480/5159/2260/5156/2321/2263/2261/1969/3791/2322/5599/5601/5602", "/")[1]
ptk_select_3 = strsplit("3717/3643/1436/3480/5159/2260/5156/2321/2263/2261/1969/3791/2322/3718", "/")[1]

ptk_select_1 = as.data.frame(ptk_select_1)
ptk_select_2 = as.data.frame(ptk_select_2)
ptk_select_3 = as.data.frame(ptk_select_3)

colnames(ptk_select_1) = colnames(ptk_select_2) = colnames(ptk_select_3) = "select"

ptk_select_1 = merge(ptk_select_1,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_2 = merge(ptk_select_2,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]
ptk_select_3 = merge(ptk_select_3,gene2_ptk, by.x = "select", by.y = "ENTREZID")[3]

stk_select_1 = strsplit("3551/1147/9252/5894/5582/5579/369/673/5578/5606/5608/7867/8491/51701", "/")[1]

stk_select_1 = as.data.frame(stk_select_1)

colnames(stk_select_1) = "select"

stk_select_1 = merge(stk_select_1,gene2_stk, by.x = "select", by.y = "ENTREZID")[3]
