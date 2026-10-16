group <- read.table("00.rawdata/group.txt",sep="\t",header = T,row.names = 1,quote = "")

analysis.list <- c("WTKO", "coHouse")
analysis.list.info <- list(
  compaired = list(
    "WTKO"=list(c("WT","KO")),
    "coHouse"=list(c("WT","KO"),c("coWT","coKO"),c("KO","coKO"),c("WT","coWT"))
  ),
  levels = list(
    "WTKO"=c("WT","KO"),
    "coHouse"=c("WT","KO","coWT","coKO")
  ),
  colors = list(
    "WTKO"=c("WT"=pal_jama()(4)[1],"KO"=pal_jama()(4)[2]),
    "coHouse"=c("WT"=pal_jama()(4)[1],"KO"=pal_jama()(4)[2],"coWT"=pal_jama()(4)[3],"coKO"=pal_jama()(4)[4])
  )
)

############################# fecal metabolites ######
dataPOS <- read.table("00.rawdata/SK_F.vs.ST_F/SK_F.vs.ST_F_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataPOS <- dataPOS[,grep("pos_",colnames(dataPOS))]
colnames(dataPOS) <- gsub("pos_STF_","sWT",colnames(dataPOS))
colnames(dataPOS) <- gsub("pos_SKF_","sKO",colnames(dataPOS))
dataNEG <- read.table("00.rawdata/SK_F.vs.ST_F/SK_F.vs.ST_F_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataNEG <- dataNEG[,grep("neg_",colnames(dataNEG))]
colnames(dataNEG) <- gsub("neg_STF_","sWT",colnames(dataNEG))
colnames(dataNEG) <- gsub("neg_SKF_","sKO",colnames(dataNEG))
######### PLS-DA ######### 
library(ropls)
library(ggplot2)
library(ggforce)
df1_plsda <- opls(t(dataPOS), group$group, orthoI = 0)
data <- as.data.frame(df1_plsda@scoreMN)
data$group = factor(c(rep("KO",10),rep("WT",10)),levels = c("WT","KO"))
data$samples = rownames(data)
x_lab <- df1_plsda@modelDF[1, "R2X"] * 100
y_lab <- df1_plsda@modelDF[2, "R2X"] * 100
col=c("#374E55","#DF8F44")
ggplot(data,aes(x=p1,y=p2,color=group,shape=group))+
  theme_bw()+
  geom_point(size=4)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("P2 (",y_lab,"%)"))+
  stat_ellipse(data=data,
               geom = "polygon",level = 0.95,
               linetype = 1,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col) +
  scale_fill_manual(values = c("#374E55","#DF8F44"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
df1_plsda <- opls(t(dataNEG), group$group, orthoI = 0)
data <- as.data.frame(df1_plsda@scoreMN)
data$group = factor(c(rep("KO",10),rep("WT",10)),levels = c("WT","KO"))
data$samples = rownames(data)
x_lab <- df1_plsda@modelDF[1, "R2X"] * 100
y_lab <- df1_plsda@modelDF[2, "R2X"] * 100
col=c("#374E55","#DF8F44")
ggplot(data,aes(x=p1,y=p2,color=group,shape=group))+
  theme_bw()+
  geom_point(size=4)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("P2 (",y_lab,"%)"))+
  stat_ellipse(data=data,
               geom = "polygon",level = 0.95,
               linetype = 1,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col) +
  scale_fill_manual(values = c("#374E55","#DF8F44"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
######### vovalno plot #########
oriPOS <- read.table("00.rawdata/SK_F.vs.ST_F/SK_F.vs.ST_F_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/SK_F.vs.ST_F/SK_F.vs.ST_F_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
p1 <- oriNEG %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("neg")
p2 <- oriPOS %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("pos")
ggarrange(p1,p2)
######### barplot #####
oriPOS <- oriPOS %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='pos',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
oriNEG <- oriNEG %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='neg',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
fecaloriSTAT <- oriSTAT
## for barplot
selClass <- (table(oriSTAT$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq) %>% dplyr::filter(Var1 != '' & Freq > 2))$Var1 %>% as.vector()
barplotIN <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
p1 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Up') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of increased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='firebrick3',low='grey')
p2 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Down') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of decreased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='navy',low='grey')
ggarrange(p1,p2,nrow = 1,ncol = 2)
######### heatmap #####
heatmapOrder <- NULL
for(midClass in selClass){
  heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% midClass)
  heat.in <- heat.in[,grep("^(SKF|STF)",colnames(heat.in))]
  mid <- pheatmap::pheatmap(heat.in,scale = "row")
  heatmapOrder <- c(heatmapOrder,mid$tree_row$labels[mid$tree_row$order])
}
heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
heat.in <- heat.in[,grep("^(SKF|STF)",colnames(heat.in))]
heat.in <- t(scale(t(heat.in),center=T,scale=T))
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
row_ha <- rowAnnotation(Class = barplotIN[heatmapOrder,"Class_II"],Type=barplotIN[heatmapOrder,"Type"]) #random for color
p1 <- Heatmap(heat.in[heatmapOrder,grep("STF_",colnames(heat.in))],  column_title = "WT", col = col_fun,cluster_rows = F,show_column_names = T,left_annotation = row_ha)
p2 <- Heatmap(heat.in[heatmapOrder,grep("SKF_",colnames(heat.in))],  column_title = "KO", col = col_fun,cluster_rows = F,show_column_names = T)
p1+p2
############################# blood metabolites ######
dataPOS <- read.table("00.rawdata/SK_B.vs.ST_B/SK_B.vs.ST_B_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataPOS <- dataPOS[,grep("pos_",colnames(dataPOS))]
colnames(dataPOS) <- gsub("pos_STB_","sWT",colnames(dataPOS))
colnames(dataPOS) <- gsub("pos_SKB_","sKO",colnames(dataPOS))
dataNEG <- read.table("00.rawdata/SK_B.vs.ST_B/SK_B.vs.ST_B_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataNEG <- dataNEG[,grep("neg_",colnames(dataNEG))]
colnames(dataNEG) <- gsub("neg_STB_","sWT",colnames(dataNEG))
colnames(dataNEG) <- gsub("neg_SKB_","sKO",colnames(dataNEG))
######### PLS-DA ######### 
# PLS-DA
library(ropls)
library(ggplot2)
library(ggforce)
df1_plsda <- opls(t(dataPOS), group$group, orthoI = 0)
data <- as.data.frame(df1_plsda@scoreMN)
data$group = factor(c(rep("KO",10),rep("WT",10)),levels = c("WT","KO"))
data$samples = rownames(data)
x_lab <- df1_plsda@modelDF[1, "R2X"] * 100
y_lab <- df1_plsda@modelDF[2, "R2X"] * 100
col=c("#374E55","#DF8F44")
ggplot(data,aes(x=p1,y=p2,color=group,shape=group))+
  theme_bw()+
  geom_point(size=4)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("P2 (",y_lab,"%)"))+
  stat_ellipse(data=data,
               geom = "polygon",level = 0.95,
               linetype = 1,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col) +
  scale_fill_manual(values = c("#374E55","#DF8F44"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
df1_plsda <- opls(t(dataNEG), group$group, orthoI = 0)
data <- as.data.frame(df1_plsda@scoreMN)
data$group = factor(c(rep("KO",10),rep("WT",10)),levels = c("WT","KO"))
data$samples = rownames(data)
x_lab <- df1_plsda@modelDF[1, "R2X"] * 100
y_lab <- df1_plsda@modelDF[2, "R2X"] * 100
col=c("#374E55","#DF8F44")
ggplot(data,aes(x=p1,y=p2,color=group,shape=group))+
  theme_bw()+
  geom_point(size=4)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color="grey")+
  geom_hline(yintercept = 0,lty="dashed",color="grey")+
  labs(x=paste0("P1 (",x_lab,"%)"),
       y=paste0("P2 (",y_lab,"%)"))+
  stat_ellipse(data=data,
               geom = "polygon",level = 0.95,
               linetype = 1,size=0.5,
               aes(fill=group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = col) +
  scale_fill_manual(values = c("#374E55","#DF8F44"))+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
######### vovalno plot #########
oriPOS <- read.table("00.rawdata/SK_B.vs.ST_B/SK_B.vs.ST_B_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/SK_B.vs.ST_B/SK_B.vs.ST_B_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
p1 <- oriNEG %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("neg")
p2 <- oriPOS %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("pos")
ggarrange(p1,p2)
######### barplot #####
oriPOS <- oriPOS %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='pos',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
oriNEG <- oriNEG %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='neg',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
bloodoriSTAT <- oriSTAT
selClass <- (table(oriSTAT$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq) %>% dplyr::filter(Var1 != ''))$Var1 %>% as.vector()
barplotIN <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
p1 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Up') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of increased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='firebrick3',low='grey')
p2 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Down') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of decreased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='navy',low='grey')
ggarrange(p1,p2,nrow = 1,ncol = 2)
######### heatmap #####
heatmapOrder <- NULL
for(midClass in selClass){
  heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% midClass)
  if(nrow(heat.in) >= 3){
    heat.in <- heat.in[,grep("^(SKB|STB)",colnames(heat.in))]
    mid <- pheatmap::pheatmap(heat.in,scale = "row")
    heatmapOrder <- c(heatmapOrder,mid$tree_row$labels[mid$tree_row$order])
  }else{
    heatmapOrder <- c(heatmapOrder,rownames(heat.in))
  }
}
heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
heat.in <- heat.in[,grep("^(SKB|STB)",colnames(heat.in))]
heat.in <- t(scale(t(heat.in),center=T,scale=T))
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
row_ha <- rowAnnotation(Class = barplotIN[heatmapOrder,"Class_II"],Type=barplotIN[heatmapOrder,"Type"]) #random for color
p1 <- Heatmap(heat.in[heatmapOrder,grep("STB_",colnames(heat.in))],  column_title = "WT", col = col_fun,cluster_rows = F,show_column_names = T,left_annotation = row_ha)
p2 <- Heatmap(heat.in[heatmapOrder,grep("SKB_",colnames(heat.in))],  column_title = "KO", col = col_fun,cluster_rows = F,show_column_names = T)
p1+p2

############################# correlation between species and metabolites ######
rownames(bloodoriSTAT) <- paste0("B_",rownames(bloodoriSTAT))
rownames(fecaloriSTAT) <- paste0("F_",rownames(fecaloriSTAT))
colnames(bloodoriSTAT) <- gsub("SKB_","SK",colnames(bloodoriSTAT))
colnames(bloodoriSTAT) <- gsub("STB_","ST",colnames(bloodoriSTAT))
colnames(fecaloriSTAT) <- gsub("SKF_","SK",colnames(fecaloriSTAT))
colnames(fecaloriSTAT) <- gsub("STF_","ST",colnames(fecaloriSTAT))
totalMBs <- rbind(bloodoriSTAT,fecaloriSTAT) %>% dplyr::filter(Class_II != '' & Class_II != "None")
totalMBs <- totalMBs[,colnames(speciesoriSTAT)]
library(dplyr)
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}

CTR.mat <- corr.test(t(totalMBs),t(speciesoriSTAT),method="spearman",adjust = "fdr")
edgeOUT <- data.frame(from=NA,to=NA,r=NA,p=NA)
i <- 1
for (mid.tax in colnames(CTR.mat$r)) {
  for (mid.rna in rownames(CTR.mat$r)) {
    if(CTR.mat$p[mid.rna,mid.tax] < 0.05 & !is.na(CTR.mat$p[mid.rna,mid.tax]) & !is.na(CTR.mat$r[mid.rna,mid.tax])){
      edgeOUT[i,] <- c(mid.rna,mid.tax,CTR.mat$r[mid.rna,mid.tax],CTR.mat$p[mid.rna,mid.tax])
      i <- i+1
    }
  }
}
edgeOUT <- edgeOUT %>% mutate(type = ifelse(r > 0,"up","down")) %>% mutate(rr=abs(as.numeric(r))) %>% dplyr::filter(rr > 0.2)
midMAT.p <- CTR.mat$p
midMAT.r <- CTR.mat$r
midMAT.p[midMAT.p > 0.05 | abs(midMAT.r) < 0.2] <- NA 
midMAT.p <- removeColsAllNa(midMAT.p) 
midMAT.p <- removeRowsAllNa(midMAT.p) 
midMAT.p[midMAT.p < 0.01] <- '+'
midMAT.p[midMAT.p < 0.05 & midMAT.p != '+'] <- '*'
midMAT.p[is.na(midMAT.p)] <- ''
heatIN <- midMAT.r[rownames(midMAT.p),colnames(midMAT.p)]
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
selMBs <- as.vector((table(edgeOUT$from) %>% as.data.frame() %>% dplyr::filter(Freq > 2))$Var1)
midAnno <- rbind(bloodoriSTAT,fecaloriSTAT)[selMBs,] %>% dplyr::filter(Class_II != '' & Class_II != "None") %>% dplyr::select(Class_II) %>% dplyr::mutate(Class_II=ifelse(Class_II == "","Others",Class_II))
pheatmap::pheatmap(t(heatIN[selMBs,]),scale = "none",cluster_row = T, cluster_col = T, border=NA, display_numbers = t(midMAT.p[selMBs,]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
orderS <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[1,1])) %>% rownames(),]),scale = "none",cluster_row = T, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[1,1])) %>% rownames(),]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p1.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p1 <- orderS$tree_col$labels[orderS$tree_col$order]
orderS <- orderS$tree_row$labels[orderS$tree_row$order]
p2 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[2,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[2,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p2.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p3 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[3,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[3,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p3.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p4 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[4,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[4,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p4.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p5 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[5,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[5,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p5.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p6 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[6,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[6,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p6.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p7 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[7,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[7,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p7.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p8 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[8,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[8,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p8.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p9 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[9,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[9,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p9.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p10 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[10,1])) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II == as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq))[10,1])) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p10.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
p11 <- pheatmap::pheatmap(t(heatIN[midAnno %>% dplyr::filter( Class_II %in% as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq) %>% dplyr::filter(Freq == 1))$Var1)) %>% rownames(),orderS]),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[midAnno %>% dplyr::filter( Class_II %in% as.vector((table(midAnno$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq) %>% dplyr::filter(Freq == 1))$Var1)) %>% rownames(),orderS]),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("species2MBs.heatmap.p11.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))


############################# correlation between species, genera, metabolites and ILs #########
ILs <- read.table("00.rawdata/COHOUSED-luminex.txt",sep = "\t",header = T,row.names = 1)
selSamples <- intersect(rownames(ILs),colnames(speciesoriSTAT))
CTR.mat <- corr.test(ILs[selSamples,],t(speciesoriSTAT[,selSamples]),method="spearman",adjust = "fdr")
edgeOUT <- data.frame(from=NA,to=NA,r=NA,p=NA)
i <- 1
for (mid.tax in colnames(CTR.mat$r)) {
  for (mid.rna in rownames(CTR.mat$r)) {
    if(CTR.mat$p[mid.rna,mid.tax] < 0.05 & !is.na(CTR.mat$p[mid.rna,mid.tax]) & !is.na(CTR.mat$r[mid.rna,mid.tax])){
      edgeOUT[i,] <- c(mid.rna,mid.tax,CTR.mat$r[mid.rna,mid.tax],CTR.mat$p[mid.rna,mid.tax])
      i <- i+1
    }
  }
}
edgeOUT <- edgeOUT %>% mutate(type = ifelse(r > 0,"up","down")) %>% mutate(rr=abs(as.numeric(r))) %>% dplyr::filter(rr > 0.2)

ILs <- read.table("00.rawdata/COHOUSED-luminex.txt",sep = "\t",header = T,row.names = 1)
selSamples <- intersect(rownames(ILs),colnames(genusoriSTAT))
CTR.mat <- corr.test(ILs[selSamples,],t(genusoriSTAT[,selSamples]),method="spearman",adjust = "fdr")
midMAT.p <- CTR.mat$p
midMAT.r <- CTR.mat$r
midMAT.p[midMAT.p > 0.05 | abs(midMAT.r) < 0.2] <- NA 
midMAT.p <- removeColsAllNa(midMAT.p) 
midMAT.p <- removeRowsAllNa(midMAT.p) 
midMAT.p[midMAT.p < 0.01] <- '+'
midMAT.p[midMAT.p < 0.05 & midMAT.p != '+'] <- '*'
midMAT.p[is.na(midMAT.p)] <- ''
heatIN <- midMAT.r[rownames(midMAT.p),colnames(midMAT.p)]
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
pheatmap::pheatmap(t(heatIN),scale = "none",cluster_row = T, cluster_col = T, border=NA, display_numbers = t(midMAT.p),fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 7, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("genus2ILs.heatmap.R0.2.pdf"),angle_col = "45",legend_breaks=c(-1,-0.5,0,0.5,1))

totalMBs <- rbind(bloodoriSTAT,fecaloriSTAT)[c(p1,p2$tree_col$labels[p2$tree_col$order],p3$tree_col$labels[p3$tree_col$order],p4$tree_col$labels[p4$tree_col$order],p5$tree_col$labels[p5$tree_col$order],p6$tree_col$labels[p6$tree_col$order],p7$tree_col$labels[p7$tree_col$order],p8$tree_col$labels[p8$tree_col$order],p9$tree_col$labels[p9$tree_col$order],p10$tree_col$labels[p10$tree_col$order],p11$tree_col$labels[p11$tree_col$order]),]
ILs <- read.table("00.rawdata/COHOUSED-luminex.txt",sep = "\t",header = T,row.names = 1)
selSamples <- intersect(rownames(ILs),colnames(totalMBs))
CTR.mat <- corr.test(ILs[selSamples,],t(totalMBs[,selSamples]),method="spearman",adjust = "fdr")
midMAT.p <- CTR.mat$p
midMAT.r <- CTR.mat$r
midMAT.p[midMAT.p > 0.05 | abs(midMAT.r) < 0.2] <- NA 
midMAT.p[midMAT.p < 0.01] <- '+'
midMAT.p[midMAT.p < 0.05 & midMAT.p != '+'] <- '*'
midMAT.p[is.na(midMAT.p)] <- ''
heatIN <- midMAT.r[rownames(midMAT.p),colnames(midMAT.p)]
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
selMBs <- rownames(totalMBs)
midAnno <- rbind(bloodoriSTAT,fecaloriSTAT)[selMBs,] %>% dplyr::filter(Class_II != '' & Class_II != "None") %>% dplyr::select(Class_II) %>% dplyr::mutate(Class_II=ifelse(Class_II == "","Others",Class_II))
pheatmap::pheatmap(heatIN[,selMBs],scale = "none",cluster_row = T, cluster_col = F, border=NA, display_numbers = midMAT.p[,selMBs],fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("MB2ILs.heatmap.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
############################# correlation between genera and metabolites #########
totalMBs <- rbind(bloodoriSTAT,fecaloriSTAT)[c(p1,p2$tree_col$labels[p2$tree_col$order],p3$tree_col$labels[p3$tree_col$order],p4$tree_col$labels[p4$tree_col$order],p5$tree_col$labels[p5$tree_col$order],p6$tree_col$labels[p6$tree_col$order],p7$tree_col$labels[p7$tree_col$order],p8$tree_col$labels[p8$tree_col$order],p9$tree_col$labels[p9$tree_col$order],p10$tree_col$labels[p10$tree_col$order],p11$tree_col$labels[p11$tree_col$order]),colnames(genusoriSTAT)]
CTR.mat <- corr.test(t(genusoriSTAT),t(totalMBs),method="spearman",adjust = "fdr")
midMAT.p <- CTR.mat$p
midMAT.r <- CTR.mat$r
midMAT.p[midMAT.p > 0.05 | abs(midMAT.r) < 0.2] <- NA 
midMAT.p[midMAT.p < 0.01] <- '+'
midMAT.p[midMAT.p < 0.05 & midMAT.p != '+'] <- '*'
midMAT.p[is.na(midMAT.p)] <- ''
heatIN <- midMAT.r[rownames(midMAT.p),colnames(midMAT.p)]
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
selMBs <- rownames(totalMBs)
midAnno <- rbind(bloodoriSTAT,fecaloriSTAT)[selMBs,] %>% dplyr::filter(Class_II != '' & Class_II != "None") %>% dplyr::select(Class_II) %>% dplyr::mutate(Class_II=ifelse(Class_II == "","Others",Class_II))
pheatmap::pheatmap(heatIN[,selMBs],scale = "none",cluster_row = T, cluster_col = F, border=NA, display_numbers = midMAT.p[,selMBs],fontsize_number = 7, number_color = "black",show_colnames = T, cellwidth = 5.5, cellheight =7,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("genus2MBs.heatmap.R0.2.pdf"),angle_col = "45",annotation_col = midAnno,legend_breaks=c(-1,-0.5,0,0.5,1))
