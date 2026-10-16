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

#############################################################################
######################### Analysis after coHoused ###########################
#############################################################################
Vir1 <- read.table("03.signFeatures/Wilcoxon.Maaslin2.s.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.coKO.VS.coWT. < 0.05) %>% dplyr::filter(coWT_Median > 0 | coKO_Median > 0) %>% rownames() 
Vir2 <- read.table("03.signFeatures/Wilcoxon.Maaslin2.g.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.coKO.VS.coWT. < 0.05) %>% dplyr::filter(coWT_Median > 0 | coKO_Median > 0) %>% rownames()
##### visualization for sign species in four groups #####
data <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[
  c('Akkermansia_muciniphila',
    'Parasutterella_excrementihominis',
    'Parabacteroides_distasonis',
    'Bacteroides_nordii',
    'Faecalibaculum_rodentium',
    'Ileibacterium_valens'),group %>% dplyr::filter(coHouse %in% c("coWT","coKO","KO","WT")) %>% rownames()] %>% tibble::rownames_to_column('ID') %>% reshape2::melt(by='ID') %>% dplyr::mutate(Group=group[variable,"coHouse"]) 
data$ID <- factor(data$ID,levels=c('Akkermansia_muciniphila',
                                   'Parasutterella_excrementihominis',
                                   'Parabacteroides_distasonis',
                                   'Bacteroides_nordii',
                                   'Faecalibaculum_rodentium',
                                   'Ileibacterium_valens'))
ggplot(data,aes(Group,log10(value),fill=Group)) + geom_boxplot() + facet_wrap('ID',scales = 'free',nrow = 2) + geom_signif(comparisons = analysis.list.info$compaired$coHouse,color="black",step_increase = 0.08,map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{'NS'}},test = wilcox.test) + theme_light() + scale_fill_manual(values = analysis.list.info$colors$coHouse)
##### visualization for sign species, coKO vs coWT #####
heat.in1 <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[Vir1,group %>% dplyr::filter(coHouse %in% c("coWT","coKO")) %>% rownames()]
heat.in2 <- read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[Vir2,group %>% dplyr::filter(coHouse %in% c("coWT","coKO")) %>% rownames()]
heat.in <- rbind(heat.in1,heat.in2)
mid <- pheatmap::pheatmap(heat.in,scale = "row")
heat.in <- t(scale(t(heat.in),center=T,scale=T))
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
p1 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "coWT") %>% rownames()],  column_title = "coWT", col = col_fun,cluster_rows = F,show_column_names = F)
p2 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "coKO") %>% rownames()],  column_title = "coKO", col = col_fun,cluster_rows = F,show_column_names = F)
p1+p2
###### visualization for sign species, KO vs WT #####
heat.in <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[rownames(speciesoriSTAT),group %>% dplyr::filter(coHouse %in% c("WT","KO") ) %>% rownames()]
mid <- pheatmap::pheatmap(heat.in,scale = "row")
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
heat.in <- t(scale(t(read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[rownames(speciesoriSTAT),group %>% dplyr::filter(coHouse %in% c("WT","KO")) %>% rownames()]),center=T,scale=T))
p1 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "WT") %>% rownames()],  column_title = "WT", col = col_fun,cluster_rows = F,show_column_names = F)
p2 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "KO") %>% rownames()],  column_title = "KO", col = col_fun,cluster_rows = F,show_column_names = F)
heat.in <- t(scale(t(read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[rownames(speciesoriSTAT),group %>% dplyr::filter(coHouse %in% c("coWT","coKO")) %>% rownames()]),center=T,scale=T))
p3 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "coWT") %>% rownames()],  column_title = "coWT", col = col_fun,cluster_rows = F,show_column_names = F)
p4 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "coKO") %>% rownames()],  column_title = "coKO", col = col_fun,cluster_rows = F,show_column_names = F)
p1+p2+p3+p4
######### genus level #########
desGenera <- read.table("03.signFeatures/Wilcoxon.Maaslin2.g.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.WTvsKO < 0.05) %>% dplyr::filter(WT_Median > 0 | KO_Median > 0) %>% rownames()
### visualization for sign species, KO vs WT #####
heat.in <- read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[desGenera,group %>% dplyr::filter(coHouse %in% c("WT","KO") ) %>% rownames()]
mid <- pheatmap::pheatmap(heat.in,scale = "row")
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
heat.in <- t(scale(t(read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[desGenera,group %>% dplyr::filter(coHouse %in% c("WT","KO")) %>% rownames()]),center=T,scale=T))
p1 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "WT") %>% rownames()],  column_title = "WT", col = col_fun,cluster_rows = F,show_column_names = F)
p2 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "KO") %>% rownames()],  column_title = "KO", col = col_fun,cluster_rows = F,show_column_names = F)
heat.in <- t(scale(t(read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[desGenera,group %>% dplyr::filter(coHouse %in% c("coWT","coKO")) %>% rownames()]),center=T,scale=T))
p3 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "coWT") %>% rownames()],  column_title = "coWT", col = col_fun,cluster_rows = F,show_column_names = F)
p4 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(coHouse == "coKO") %>% rownames()],  column_title = "coKO", col = col_fun,cluster_rows = F,show_column_names = F)
p1+p2+p3+p4

######### fecal metabolites ##########
oriPOS <- read.table("00.rawdata/CK_F.vs.CT_F/CK_F.vs.CT_F_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/CK_F.vs.CT_F/CK_F.vs.CT_F_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
p1 <- oriNEG %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("neg")
p2 <- oriPOS %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("pos")
######### barplot
oriPOS <- oriPOS %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='pos',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
oriNEG <- oriNEG %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='neg',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
fecaloriSTATcoHoused <- oriSTAT
## for barplot
selClass <- (table(oriSTAT$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq) %>% dplyr::filter(Var1 != ''))$Var1 %>% as.vector()
barplotIN <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
p1 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Up') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of increased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='firebrick3',low='grey')
p2 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Down') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of decreased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='navy',low='grey')
ggsave("coHoused/fecal.signClass2.barplot.pdf",ggarrange(p1,p2,nrow = 1,ncol = 2),width = 10,height = 3)
######### heatmap 
heatmapOrder <- NULL
for(midClass in selClass){
  heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% midClass)
  if(nrow(heat.in) >= 3){
    heat.in <- heat.in[,grep("^(CKF|CTF)",colnames(heat.in))]
    mid <- pheatmap::pheatmap(heat.in,scale = "row")
    heatmapOrder <- c(heatmapOrder,mid$tree_row$labels[mid$tree_row$order])
  }else{
    heatmapOrder <- c(heatmapOrder,rownames(heat.in))
  }
}
heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% selClass) 
heat.in <- heat.in[,grep("^(CKF|CTF)",colnames(heat.in))]
heat.in <- t(scale(t(heat.in),center=T,scale=T))
#heat.in <- oriSTAT %>% dplyr::filter(Class_II == '') 
#heat.in <- heat.in[,grep("^(CKF|CTF)",colnames(heat.in))]
#mid <- pheatmap::pheatmap(heat.in,scale = "row")
#heatmapOrder <- mid$tree_row$labels[mid$tree_row$order]
#heat.in <- t(scale(t(heat.in),center=T,scale=T))
#row_ha <- rowAnnotation(Class = rep("Others",nrow(heat.in)),Type=oriSTAT[heatmapOrder,"Type"])
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
row_ha <- rowAnnotation(Class = factor(barplotIN[heatmapOrder,"Class_II"],levels = selClass),Type=barplotIN[heatmapOrder,"Type"])
p1 <- Heatmap(heat.in[heatmapOrder,grep("CTF_",colnames(heat.in))],  column_title = "coWT", col = col_fun,cluster_rows = F,show_column_names = T,left_annotation = row_ha)
p2 <- Heatmap(heat.in[heatmapOrder,grep("CKF_",colnames(heat.in))],  column_title = "coKO", col = col_fun,cluster_rows = F,show_column_names = T)
pdf("coHoused/fecal.signClass2.heatmap.pdf",width = 10,height = 8)
p1+p2
dev.off()
######### PLS-DA for fecal ######### 
dataPOS <- read.table("00.rawdata/CK_F.vs.CT_F/CK_F.vs.CT_F_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataPOS <- dataPOS[,grep("pos_",colnames(dataPOS))]
colnames(dataPOS) <- gsub("pos_CTF_","coWT",colnames(dataPOS))
colnames(dataPOS) <- gsub("pos_CKF_","coKO",colnames(dataPOS))
dataNEG <- read.table("00.rawdata/CK_F.vs.CT_F/CK_F.vs.CT_F_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataNEG <- dataNEG[,grep("neg_",colnames(dataNEG))]
colnames(dataNEG) <- gsub("neg_CTF_","coWT",colnames(dataNEG))
colnames(dataNEG) <- gsub("neg_CKF_","coKO",colnames(dataNEG))

# PLS-DA
# same as previous code
######### blood metabolites ##########
oriPOS <- read.table("00.rawdata/CK_B.vs.CT_B/CK_B.vs.CT_B_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/CK_B.vs.CT_B/CK_B.vs.CT_B_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
p1 <- oriNEG %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("neg")
p2 <- oriPOS %>% dplyr::mutate(Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign')) %>% ggplot(aes(log2FC,-log10(Pvalue),size=VIP,color=Type)) + scale_color_manual(values = c("#3B4992BF","#b2b2b2BF","#EE0000BF")) + scale_size("size_area", range = c(-5, 5)) + theme_bw() + geom_hline(yintercept = c(-log10(0.05)),linetype=3) + geom_vline(xintercept = c(1,-1),linetype=3)  + geom_point() + ggtitle("pos")
######### barplot
oriPOS <- oriPOS %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='pos',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
oriNEG <- oriNEG %>% dplyr::filter(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1) %>% dplyr::mutate(Mode='neg',Type=ifelse(VIP > 1 & Pvalue < 0.05 & abs(log2FC) > 1,ifelse(log2FC > 0,'Up','Down'),'noSign'))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
bloodoriSTATcoHoused <- oriSTAT
## for barplot
selClass <- (table(oriSTAT$Class_II) %>% as.data.frame() %>% dplyr::arrange(-Freq) %>% dplyr::filter(Var1 != ''))$Var1 %>% as.vector()
barplotIN <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
p1 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Up') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of increased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='firebrick3',low='grey')
p2 <- table(as.vector(barplotIN$Class_II),as.vector(barplotIN$Type)) %>% as.data.frame() %>% dplyr::filter(Var2 == 'Down') %>% ggplot(aes(Var1,weight=Freq,fill=Freq)) + geom_bar() + coord_flip() + theme_classic() + ylab("Number of decreased MBs") + xlim(rev(selClass)) + xlab("") + scale_fill_gradient(high='navy',low='grey')
######### heatmap 
heatmapOrder <- NULL
for(midClass in selClass){
  heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% midClass)
  if(nrow(heat.in) >= 3){
    heat.in <- heat.in[,grep("^(CKB|CTB)",colnames(heat.in))]
    mid <- pheatmap::pheatmap(heat.in,scale = "row")
    heatmapOrder <- c(heatmapOrder,mid$tree_row$labels[mid$tree_row$order])
  }else{
    heatmapOrder <- c(heatmapOrder,rownames(heat.in))
  }
}
heat.in <- oriSTAT %>% dplyr::filter(Class_II %in% selClass)
heat.in <- heat.in[,grep("^(CKB|CTB)",colnames(heat.in))]
heat.in <- t(scale(t(heat.in),center=T,scale=T))
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
row_ha <- rowAnnotation(Class = factor(barplotIN[heatmapOrder,"Class_II"],levels = selClass),Type=barplotIN[heatmapOrder,"Type"]) #random for color
p1 <- Heatmap(heat.in[heatmapOrder,grep("CTB_",colnames(heat.in))],  column_title = "coWT", col = col_fun,cluster_rows = F,show_column_names = T,left_annotation = row_ha)
p2 <- Heatmap(heat.in[heatmapOrder,grep("CKB_",colnames(heat.in))],  column_title = "coKO", col = col_fun,cluster_rows = F,show_column_names = T)
p1+p2

######### PLS-DA for blood ######### 
dataPOS <- read.table("00.rawdata/CK_B.vs.CT_B/CK_B.vs.CT_B_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataPOS <- dataPOS[,grep("pos_",colnames(dataPOS))]
colnames(dataPOS) <- gsub("pos_CTB_","coWT",colnames(dataPOS))
colnames(dataPOS) <- gsub("pos_CKB_","coKO",colnames(dataPOS))
dataNEG <- read.table("00.rawdata/CK_B.vs.CT_B/CK_B.vs.CT_B_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
dataNEG <- dataNEG[,grep("neg_",colnames(dataNEG))]
colnames(dataNEG) <- gsub("neg_CTB_","coWT",colnames(dataNEG))
colnames(dataNEG) <- gsub("neg_CKB_","coKO",colnames(dataNEG))

# PLS-DA
# same as previous code
######### WTvsKO   VS. coWTvscoKO in fecal metabolites #######
### vennplot
vennPlot <- venn.diagram(list(WTvsKO=gsub("^F_","",rownames(fecaloriSTAT)),coWTvscoKO=rownames(fecaloriSTATcoHoused)),filename = NULL,fill=c("navy","firebrick3"))
grid.draw(vennPlot)
######### WTvsKO   VS. coWTvscoKO in blood metabolites #######
### vennplot
vennPlot <- venn.diagram(list(WTvsKO=gsub("^B_","",rownames(bloodoriSTAT)),coWTvscoKO=rownames(bloodoriSTATcoHoused)),filename = NULL,fill=c("navy","firebrick3"))
grid.draw(vennPlot)

############################# correlation between species and metabolites ######
speciesoriSTAT <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[read.table("03.signFeatures/Wilcoxon.Maaslin2.s.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.coKO.VS.coWT. < 0.05) %>% dplyr::filter(coWT_Median > 0 | coKO_Median > 0) %>% rownames(),group %>% dplyr::filter(coHouse %in% c("coWT","coKO")) %>% rownames()]
rownames(bloodoriSTATcoHoused) <- paste0("B_",rownames(bloodoriSTATcoHoused))
rownames(fecaloriSTATcoHoused) <- paste0("F_",rownames(fecaloriSTATcoHoused))
colnames(bloodoriSTATcoHoused) <- gsub("CKB_","CK",colnames(bloodoriSTATcoHoused))
colnames(bloodoriSTATcoHoused) <- gsub("CTB_","CT",colnames(bloodoriSTATcoHoused))
colnames(fecaloriSTATcoHoused) <- gsub("CKF_","CK",colnames(fecaloriSTATcoHoused))
colnames(fecaloriSTATcoHoused) <- gsub("CTF_","CT",colnames(fecaloriSTATcoHoused))
totalMBs <- rbind(bloodoriSTATcoHoused,fecaloriSTATcoHoused) # %>% dplyr::filter(Class_II != '' & Class_II != "None")
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
selMBs <- rownames(heatIN)
midAnno <- rbind(bloodoriSTATcoHoused,fecaloriSTATcoHoused)[selMBs,] %>% dplyr::filter(Class_II != '' & Class_II != "None") %>% dplyr::select(Class_II) %>% dplyr::mutate(Class_II=ifelse(Class_II == "","Others",Class_II))
pheatmap::pheatmap(heatIN[selMBs,],scale = "none",cluster_row = T, cluster_col = T, border=NA, display_numbers = midMAT.p[selMBs,],fontsize_number = 5, number_color = "black",show_colnames = T, cellwidth = 8, cellheight =5.8,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("coHoused/species2MBs.heatmap.R0.2.pdf"),angle_col = "45",legend_breaks=c(-1,-0.5,0,0.5,1)) #annotation_col = midAnno,

############################# correlation between genus and metabolites ######
genusoriSTAT <- read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[read.table("03.signFeatures/Wilcoxon.Maaslin2.g.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.coKO.VS.coWT. < 0.05) %>% dplyr::filter(coWT_Median > 0 | coKO_Median > 0) %>% rownames(),group %>% dplyr::filter(coHouse %in% c("coWT","coKO")) %>% rownames()]
totalMBs <- rbind(bloodoriSTATcoHoused,fecaloriSTATcoHoused) # %>% dplyr::filter(Class_II != '' & Class_II != "None")
totalMBs <- totalMBs[,colnames(genusoriSTAT)]
library(dplyr)
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}

CTR.mat <- corr.test(t(totalMBs),t(genusoriSTAT),method="spearman",adjust = "fdr")
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
selMBs <- rownames(heatIN)
pheatmap::pheatmap(heatIN[selMBs,],scale = "none",cluster_row = T, cluster_col = T, border=NA, display_numbers = midMAT.p[selMBs,],fontsize_number = 5, number_color = "black",show_colnames = T, cellwidth = 8, cellheight =5.8,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("coHoused/genus2MBs.heatmap.R0.2.pdf"),angle_col = "45",legend_breaks=c(-1,-0.5,0,0.5,1)) 

############################# correlation between 6 species and all metabolites in four groups #####
speciesoriSTAT <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[
  c('Akkermansia_muciniphila',
    'Parasutterella_excrementihominis',
    'Parabacteroides_distasonis',
    'Bacteroides_nordii',
    'Faecalibaculum_rodentium',
    'Ileibacterium_valens'),group %>% dplyr::filter(coHouse %in% c("coWT","coKO","KO","WT")) %>% rownames()]
oriPOS <- read.table("00.rawdata/SK_B.vs.ST_B/SK_B.vs.ST_B_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/SK_B.vs.ST_B/SK_B.vs.ST_B_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
bloodoriSTAT <- oriSTAT[,c(grep("SK|ST",colnames(oriSTAT)))]
oriPOS <- read.table("00.rawdata/CK_B.vs.CT_B/CK_B.vs.CT_B_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/CK_B.vs.CT_B/CK_B.vs.CT_B_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
bloodafterSTAT <- oriSTAT[,c(grep("CK|CT",colnames(oriSTAT)),10)]
bloodoriSTAT <- cbind(bloodoriSTAT,bloodafterSTAT)

oriPOS <- read.table("00.rawdata/SK_F.vs.ST_F/SK_F.vs.ST_F_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/K_F.vs.ST_F/SK_F.vs.ST_F_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
fecaloriSTAT <- oriSTAT[,c(grep("SK|ST",colnames(oriSTAT)))]
oriPOS <- read.table("00.rawdata/CK_F.vs.CT_F/CK_F.vs.CT_F_pos.xls",sep = "\t",header = T,row.names = 2,quote = "")
oriNEG <- read.table("00.rawdata/CK_F.vs.CT_F/CK_F.vs.CT_F_neg.xls",sep = "\t",header = T,row.names = 2,quote = "")
colnames(oriPOS) <- gsub("pos_","",colnames(oriPOS))
colnames(oriNEG) <- gsub("neg_","",colnames(oriNEG))
oriSTAT <- rbind(oriNEG,oriPOS)
fecalafterSTAT <- oriSTAT[,c(grep("CK|CT",colnames(oriSTAT)),10)]
fecaloriSTAT <- cbind(fecaloriSTAT,fecalafterSTAT)

rownames(bloodoriSTAT) <- paste0("B_",rownames(bloodoriSTAT))
colnames(bloodoriSTAT) <- gsub("SKB_","SK",colnames(bloodoriSTAT))
colnames(bloodoriSTAT) <- gsub("STB_","ST",colnames(bloodoriSTAT))
colnames(bloodoriSTAT) <- gsub("CKB_","CK",colnames(bloodoriSTAT))
colnames(bloodoriSTAT) <- gsub("CTB_","CT",colnames(bloodoriSTAT))

rownames(fecaloriSTAT) <- paste0("F_",rownames(fecaloriSTAT))
colnames(fecaloriSTAT) <- gsub("SKF_","SK",colnames(fecaloriSTAT))
colnames(fecaloriSTAT) <- gsub("STF_","ST",colnames(fecaloriSTAT))
colnames(fecaloriSTAT) <- gsub("CKF_","CK",colnames(fecaloriSTAT))
colnames(fecaloriSTAT) <- gsub("CTF_","CT",colnames(fecaloriSTAT))

totalMBs <- rbind(bloodoriSTAT,fecaloriSTAT) %>% as.data.frame() %>% dplyr::filter(Class_II != '' & Class_II != "None")
totalMBsNames <- totalMBs %>% dplyr::select(Class_II) %>% tibble::rownames_to_column('Name')
rownames(totalMBsNames) <- seq(1:nrow(totalMBs))
totalMBs <- totalMBs[,colnames(speciesoriSTAT)]
rownames(totalMBs) <- seq(1:nrow(totalMBs))


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
edgeOUT <- edgeOUT %>% mutate(type = ifelse(r > 0,"up","down")) %>% mutate(rr=abs(as.numeric(r))) %>% dplyr::filter(rr > 0.4)
edgeOUT$from %>% unique() %>% length()

midMAT.p <- CTR.mat$p
midMAT.r <- CTR.mat$r
midMAT.p[midMAT.p > 0.05 | abs(midMAT.r) < 0.5] <- NA 
midMAT.p <- removeColsAllNa(midMAT.p) 
midMAT.p <- removeRowsAllNa(midMAT.p) 
midMAT.p[midMAT.p < 0.01] <- '+'
midMAT.p[midMAT.p < 0.05 & midMAT.p != '+'] <- '*'
midMAT.p[is.na(midMAT.p)] <- ''
heatIN <- midMAT.r[rownames(midMAT.p),colnames(midMAT.p)]
bk <- c(seq(-1,0,by=0.01),seq(0.01,1,by=0.01))
selMBs <- rownames(heatIN)
rownames(heatIN) <- totalMBsNames[rownames(heatIN),'Name']
midAnno <- totalMBsNames %>% tibble::column_to_rownames('Name') %>% dplyr::select(Class_II) 
pheatmap::pheatmap(t(heatIN),scale = "none",cluster_row = F, cluster_col = T, border=NA, display_numbers = t(midMAT.p[selMBs,]),fontsize_number = 5, number_color = "black",show_colnames = T, cellwidth = 8, cellheight =5.8,fontsize = 5,color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),breaks=bk,filename = paste0("6species2MBs.heatmap.R0.5.pdf"),angle_col = "45",legend_breaks=c(-1,-0.5,0,0.5,1))#,annotation_col = midAnno)
