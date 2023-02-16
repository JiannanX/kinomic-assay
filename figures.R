# figure 2. heatmap for FC + P
getwd()
setwd("/Users/nicole/Desktop/ic/results/results1")
options(encoding = "UTF-8")

library(tidyr)
library(dplyr)
library(circlize)
library(grid)
library(pheatmap)

#figure for the heatmap of signal intensity 
# READ CSV 
d1=read.csv("FC4.csv", header = T)
d2=read.csv("P4.csv",header = T)

d1<-na.omit(d1)
d2<-na.omit(d2)

rownames(d1)=d1$ID
rownames(d2)=d2$ID

d1=d1[,c(-1)]
d2=d2[,c(-1)]

star2 = d2
for (i in 1:121){for(j in 1:3){if(d2[i,j]<0.05){star2[i,j]<-"•";}else{star2[i,j]<-"";}}}
for (i in 1:121){for(j in 1:3){if(d2[i,j]<0.01){star2[i,j]<-"••";}}}
for (i in 1:121){for(j in 1:3){if(d2[i,j]<0.001){star2[i,j]<-"•••";}}}

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
tiff("stk2.tiff",width = 4000, height = 8000, res = 400) #ptk1: height=13000, stk1:9000, ptk2:8000
pheatmap(d1,scale = "row", 
         
         color = c(colorRampPalette(colors = c("yellow","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","green"))(length(bk)/2)),
         legend_breaks=seq(-2,2,2),
         breaks=bk,
         angle_col = 45, fontsize_row = 9, fontsize_col = 12, 
         cellwidth = 100, cellheight = 9,
         
         treeheight_col = 10, 
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         cluster_row = F,
         cluster_cols = F,
         show_colnames = F,
         show_rownames = T,
         
         display_numbers = star2,
fontsize_number = 26, number_color = "black")

dev.off()


# FIGURE 3 heatmap for ex4，phe & asp
setwd("/Users/nicole/Desktop/11.6/ptk")
library(tidyr)
library(dplyr)
library(circlize)
library(grid)
library(pheatmap)

database=read.csv("ptkdataset1.csv", header = T) #database
database = database[-which(database$rank == "Kinase 26:"),]
database = database[-which(database$rank == "Kinase 27:"),]
database = database[-which(database$rank == "Kinase 28:"),]
database = database[-which(database$rank == "Kinase 29:"),]
database = database[-which(database$rank == "Kinase 30:"),]
database = database[-which(database$rank == "Kinase 31:"),]
database = database[-which(database$rank == "Kinase 32:"),]
database = database[-which(database$rank == "Kinase 33:"),]
database = database[-which(database$rank == "Kinase 34:"),]
database = database[-which(database$rank == "Kinase 35:"),]
database = database[-which(database$rank == "Kinase 36:"),]
database = database[-which(database$rank == "Kinase 37:"),]
database = database[-which(database$rank == "Kinase 38:"),]
database = database[-which(database$rank == "Kinase 39:"),]
database = database[-which(database$rank == "Kinase 40:"),]
database = database[-which(database$rank == "Kinase 41:"),]
database = database[-which(database$rank == "Kinase 42:"),]
database = database[-which(database$rank == "Kinase 43:"),]
database = database[-which(database$rank == "Kinase 44:"),]
database = database[-which(database$rank == "Kinase 45:"),]
database = database[-which(database$rank == "Kinase 46:"),]
database = database[-which(database$rank == "Kinase 47:"),]
database = database[-which(database$rank == "Kinase 48:"),]
database = database[-which(database$rank == "Kinase 49:"),]
database = database[-which(database$rank == "Kinase 50:"),]



d1=read.csv("endoc101.csv", header = T)
d2=read.csv("endoc102.csv", header = T)
d3=read.csv("endoc103.csv", header = T)
d1<-na.omit(d1)
d2<-na.omit(d2)
d3<-na.omit(d3)


d1<- d1%>% as_tibble() %>% 
  separate_rows(sites, sep = ";")

d2<- d2%>% as_tibble() %>% 
  separate_rows(sites, sep = ";")

d3<- d3%>% as_tibble() %>% 
  separate_rows(sites, sep = ";")

px1 = merge(database, d1)
px2 = merge(database, d2)
px3 = merge(database, d3)

px1 = px1[order(px1$kinase),]
px2 = px2[order(px2$kinase),]
px3 = px3[order(px3$kinase),]

ex = px1[,c(6,1,2,3,8,9)]
phe = px2[,c(6,1,2,3,8,9)]
asp = px3[,c(6,1,2,3,8,9)]
ex<-na.omit(ex)
phe<-na.omit(phe)
asp<-na.omit(asp)

colnames(ex)=c("kinase","sub","sites","peptides","score","logFC")
colnames(phe)=c("kinase","sub","sites","peptides","score","logFC")
colnames(asp)=c("kinase","sub","sites","peptides","score","logFC")

ptkex<-group_by(ex,kinase)
ptkex<-summarise(ptkex,Sw=weighted.mean(logFC,score))

ptkphe<-group_by(phe,kinase)
ptkphe<-summarise(ptkphe,Sw=weighted.mean(logFC,score))

ptkasp<-group_by(asp,kinase)
ptkasp<-summarise(ptkasp,Sw=weighted.mean(logFC,score))


kinase.list1 = as.vector(ex$kinase)
kinase.list1 = as.matrix(table(kinase.list1))

kinase.list2 = as.vector(phe$kinase)
kinase.list2 = as.matrix(table(kinase.list2))

kinase.list3 = as.vector(asp$kinase)
kinase.list3 = as.matrix(table(kinase.list3))


#calculate z
#ex
score1 = ptkex
score1$Enrichment = score1$Sw/abs(mean(d1$ex4, na.rm=T))
score1$m = kinase.list1
score1$z = ((score1$Sw- mean(d1$ex4, na.rm=T))*sqrt(score1$m))/sd(d1$ex4, na.rm=T)
score1$p.value = pnorm(-abs(score1$z)) # 1-tailed p-value
score1$FDR = p.adjust(score1$p.value, method="fdr")
#phe
score2 = ptkphe
score2$Enrichment = score2$Sw/abs(mean(d2$phe, na.rm=T))
score2$m = kinase.list2
score2$z = ((score2$Sw- mean(d2$phe, na.rm=T))*sqrt(score2$m))/sd(d2$phe, na.rm=T)
score2$p.value = pnorm(-abs(score2$z)) # 1-tailed p-value
score2$FDR = p.adjust(score2$p.value, method="fdr")
#asp
score3 = ptkasp
score3$Enrichment = score3$Sw/abs(mean(d3$asp, na.rm=T))
score3$m = kinase.list3
score3$z = ((score3$Sw- mean(d3$asp, na.rm=T))*sqrt(score3$m))/sd(d3$asp, na.rm=T)
score3$p.value = pnorm(-abs(score3$z)) # 1-tailed p-value
score3$FDR = p.adjust(score3$p.value, method="fdr")

#filtered
score.123 = union(score1$kinase[score1$p.value < 0.05],score2$kinase[score2$p.value < 0.05])
score.123 = union(score.123,score3$kinase[score3$p.value < 0.05])
score.all.z = cbind(score1[score1$kinase %in% score.123,c(1,5)], score2[score2$kinase %in% score.123,5], score3[score3$kinase %in% score.123,5])
colnames(score.all.z) = c("kinase","ex4","phe1","asp3")
score.all.p = cbind(score1[score1$kinase %in% score.123,c(1,6)], score2[score2$kinase %in% score.123,6], score3[score3$kinase %in% score.123,6])
colnames(score.all.p) = c("kinase","ex4","phe1","asp3")

write.csv(score.all.z,"stk10all.csv")

rownames(score.all.z)=score.all.z$kinase
rownames(score.all.p)=score.all.p$kinase
z1=score.all.z[,c(-1)]
new=score.all.p[,c(-1)]


star2 = new
for (i in 1:113){for(j in 1:3){if(new[i,j]<0.05){star2[i,j]<-"•";}else{star2[i,j]<-"";}}}
for (i in 1:113){for(j in 1:3){if(new[i,j]<0.01){star2[i,j]<-"••";}}}
for (i in 1:113){for(j in 1:3){if(new[i,j]<0.001){star2[i,j]<-"•••";}}}


bk <- c(seq(-6,-0.1,by=0.01),seq(0,6,by=0.01))

tiff("STK10.tiff",width = 1600, height = 6000, res = 400)
pheatmap(z1,scale = "none",
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-6,6,2),
         breaks=bk,
         angle_col = 0, fontsize_row = 8, fontsize_col = 12, 
         cellwidth = 30, cellheight = 8,
        
         treeheight_col = 10, 
         cluster_cols = F,
         cluster_row = FALSE, display_numbers = star2,
         fontsize_number = 10 , number_color = "black")
dev.off()






# phe vs asp

setwd("/Users/nicole/Desktop/11.6/phe&asp")
library(tidyr)
library(dplyr)
library(circlize)
library(grid)
library(pheatmap)

database=read.csv("stkdataset1.csv", header = T) #database
database = database[-which(database$rank == "Kinase 26:"),]
database = database[-which(database$rank == "Kinase 27:"),]
database = database[-which(database$rank == "Kinase 28:"),]
database = database[-which(database$rank == "Kinase 29:"),]
database = database[-which(database$rank == "Kinase 30:"),]
database = database[-which(database$rank == "Kinase 31:"),]
database = database[-which(database$rank == "Kinase 32:"),]
database = database[-which(database$rank == "Kinase 33:"),]
database = database[-which(database$rank == "Kinase 34:"),]
database = database[-which(database$rank == "Kinase 35:"),]
database = database[-which(database$rank == "Kinase 36:"),]
database = database[-which(database$rank == "Kinase 37:"),]
database = database[-which(database$rank == "Kinase 38:"),]
database = database[-which(database$rank == "Kinase 39:"),]
database = database[-which(database$rank == "Kinase 40:"),]
database = database[-which(database$rank == "Kinase 41:"),]
database = database[-which(database$rank == "Kinase 42:"),]
database = database[-which(database$rank == "Kinase 43:"),]
database = database[-which(database$rank == "Kinase 44:"),]
database = database[-which(database$rank == "Kinase 45:"),]
database = database[-which(database$rank == "Kinase 46:"),]
database = database[-which(database$rank == "Kinase 47:"),]
database = database[-which(database$rank == "Kinase 48:"),]
database = database[-which(database$rank == "Kinase 49:"),]
database = database[-which(database$rank == "Kinase 50:"),]



d1=read.csv("s1phe.csv", header = T)
d2=read.csv("s1asp.csv", header = T)

d1<-na.omit(d1)
d2<-na.omit(d2)


d1<- d1%>% as_tibble() %>% 
  separate_rows(sites, sep = ";")

d2<- d2%>% as_tibble() %>% 
  separate_rows(sites, sep = ";")


px1 = merge(database, d1)
px2 = merge(database, d2)


px1 = px1[order(px1$kinase),]
px2 = px2[order(px2$kinase),]

phe = px1[,c(6,1,2,3,8,9)]
asp = px2[,c(6,1,2,3,8,9)]

phe<-na.omit(phe)
asp<-na.omit(asp)


colnames(phe)=c("kinase","sub","sites","peptides","score","logFC")
colnames(asp)=c("kinase","sub","sites","peptides","score","logFC")


ptkphe<-group_by(phe,kinase)
ptkphe<-summarise(ptkphe,Sw=weighted.mean(logFC,score))

ptkasp<-group_by(asp,kinase)
ptkasp<-summarise(ptkasp,Sw=weighted.mean(logFC,score))


kinase.list2 = as.vector(phe$kinase)
kinase.list2 = as.matrix(table(kinase.list2))

kinase.list3 = as.vector(asp$kinase)
kinase.list3 = as.matrix(table(kinase.list3))


#calculate z

#phe
score2 = ptkphe
score2$Enrichment = score2$Sw/abs(mean(d1$phe, na.rm=T))
score2$m = kinase.list2
score2$z = ((score2$Sw- mean(d1$phe, na.rm=T))*sqrt(score2$m))/sd(d1$phe, na.rm=T)
score2$p.value = pnorm(-abs(score2$z)) # 1-tailed p-value
score2$FDR = p.adjust(score2$p.value, method="fdr")
#asp
score3 = ptkasp
score3$Enrichment = score3$Sw/abs(mean(d2$asp, na.rm=T))
score3$m = kinase.list3
score3$z = ((score3$Sw- mean(d2$asp, na.rm=T))*sqrt(score3$m))/sd(d2$asp, na.rm=T)
score3$p.value = pnorm(-abs(score3$z)) # 1-tailed p-value
score3$FDR = p.adjust(score3$p.value, method="fdr")

#filtered
score.23 = union(score2$kinase[score2$p.value < 0.05],score3$kinase[score3$p.value < 0.05])
score.all.z = cbind(score2[score2$kinase %in% score.23,c(1,5)], score3[score3$kinase %in% score.23,5])
colnames(score.all.z) = c("kinase","phe1","asp3")
score.all.p = cbind(score2[score2$kinase %in% score.23,c(1,6)], score3[score3$kinase %in% score.23,6])
colnames(score.all.p) = c("kinase","phe1","asp3")

write.csv(score.all.z,"ptk10all.csv")

rownames(score.all.z)=score.all.z$kinase
rownames(score.all.p)=score.all.p$kinase
z1=score.all.z[,c(-1)]
new=score.all.p[,c(-1)]

star2 = new
for (i in 1:50){for(j in 1:2){if(new[i,j]<0.05){star2[i,j]<-"•";}else{star2[i,j]<-"";}}}
for (i in 1:50){for(j in 1:2){if(new[i,j]<0.01){star2[i,j]<-"••";}}}
for (i in 1:50){for(j in 1:2){if(new[i,j]<0.001){star2[i,j]<-"•••";}}}


bk <- c(seq(-6,-0.1,by=0.01),seq(0,6,by=0.01))
tiff("STK10.tiff",width = 1500, height = 2600, res = 400)
pheatmap(z1,scale = "none",
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-6,6,2),
         breaks=bk,
         angle_col = 0, fontsize_row = 8, fontsize_col = 12, 
         cellwidth = 30, cellheight = 8,
         
         treeheight_col = 10, 
         cluster_cols = F,
         cluster_row = FALSE, display_numbers = star2,
         fontsize_number = 10 , number_color = "black")
dev.off()
