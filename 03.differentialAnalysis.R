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

############################## FUNCTIONS ###############################
group_mean <- function(group.subgroup,data,colInfo="Group"){
  tmpg=unique(group.subgroup[,colInfo])
  tmpn=length(tmpg)
  rows=nrow(data)
  result <- matrix(nrow=rows,ncol=tmpn*3)
  name <- NULL
  for(i in tmpg){ name <- c(name,paste(i,'Median',sep='_'),paste(i,'Mean',sep='_'),paste(i,'SD',sep='_'))    }
  colnames(result)=name
  for(i in tmpg){
    tmpdata=data[,rownames(group.subgroup[which(group.subgroup[,colInfo]==i),,drop=F])]
    name <- c(paste(i,'Median',sep='_'),paste(i,'Mean',sep='_'),paste(i,'SD',sep='_'))
    for (j in 1:rows){
      tmp.median=median(as.numeric(tmpdata[j,]))
      tmp.mean=mean(as.numeric(tmpdata[j,]))
      tmp.SD=sd(as.numeric(tmpdata[j,]))
      result[j,name] <- c(tmp.median,tmp.mean,tmp.SD)
    }
  }
  rownames(result) <- rownames(data)
  result <- data.frame(result)
  return(result)
}
wilcoxon.FDR <- function(indir,outdir.main){
  if(! file.exists(outdir.main)){dir.create(outdir.main)}
  
  for (ana.name in analysis.list) {
    class=c('p','g','s')
    for (cls in class){
      
      if(! file.exists(paste0(outdir.main,ana.name,"/"))){dir.create(paste0(outdir.main,ana.name,"/"))}
      
      outfile=paste(outdir.main,ana.name,'/Wilcoxon.',cls,'.xls',sep='')
      infile=paste(indir,'/',cls,'_tax_rel.sort.screen.txt',sep='')
      gfile='00.rawdata/group.txt'
      data <- read.table(infile,sep="\t",header=T,row.names=1)
      group <- read.table(gfile,sep="\t",header=T,row.names=1,stringsAsFactors=F)
      if(ana.name != "SampleType"){group <- group[which(group[,ana.name] != "NA"),]}
      group <- group[intersect(rownames(group),colnames(data)),]
      data=data[,rownames(group)]
      tmpg=unique(group[,ana.name])
      tmpn=length(tmpg)
      rows=nrow(data)
      #statistics
      result <- group_mean(group,data,ana.name)
      for (i in 1:(tmpn-1)) {
        for (j in (i+1):tmpn){
          #wilcoxon test
          sample1=rownames(group[which(group[,ana.name]==tmpg[i]),,drop=F])
          sample2=rownames(group[which(group[,ana.name]==tmpg[j]),,drop=F])
          wilcox.01.p <- apply(data,1,function(x) {if(sum(x[c(sample1,sample2)]) == 0){return('NA')}else{test <- wilcox.test(x[sample1],x[sample2],conf.int = T);p <- test$p.value;return(p);}}) %>% as.numeric()
          wilcox.01.low <- apply(data,1,function(x) {if(sum(x[c(sample1,sample2)]) == 0){return('NA')}else{test <- wilcox.test(x[sample1],x[sample2],conf.int = T);p <- as.vector(test$conf.int)[1];return(p);}}) %>% as.numeric()
          wilcox.01.high <- apply(data,1,function(x) {if(sum(x[c(sample1,sample2)]) == 0){return('NA')}else{test <- wilcox.test(x[sample1],x[sample2],conf.int = T);p <- as.vector(test$conf.int)[2];return(p);}}) %>% as.numeric()
          # wilcox.01.q <- p.adjust(wilcox.01.p, method = "fdr")
          wilcox.01.q <- rep("NA",length(wilcox.01.p))
          wilcox.01.q[which(wilcox.01.p < 0.05)] <- p.adjust(wilcox.01.p[which(wilcox.01.p < 0.05)], method = "fdr")
          wilname.p <- data.frame(p=wilcox.01.p,q=wilcox.01.q,low=wilcox.01.low,high=wilcox.01.high)
          colnames(wilname.p) <- paste0(c("wilcox.test.p","wilcox.test.q","low","high"),paste0("(",tmpg[i],' VS ',tmpg[j],")"))
          result <- cbind(result,wilname.p)
        }
      }
      write.table(result,file=outfile,sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
    }
  }
}
############################## diff tax ################################
wilcoxon.FDR(indir="00.rawdata/",outdir.main="03.signFeatures/")
######### species level #########
selSpecies <- read.table("03.signFeatures/Wilcoxon.s.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.WTvsKO < 0.05) %>% dplyr::filter(WT_Median > 0 | KO_Median > 0) %>% rownames()
read.table("03.signFeatures/Wilcoxon.s.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.WTvsKO < 0.05) %>% dplyr::filter(WT_Median > 0 | KO_Median > 0) %>% dplyr::mutate(log2FC=ifelse(KO_Mean > WT_Mean,log2(KO_Mean/WT_Mean),-log2(WT_Mean/KO_Mean))) %>% tibble::rownames_to_column('ID') %>%  dplyr::mutate(Col=ifelse(log2FC > 0,'red','blue')) %>% ggplot(aes(weight=log2FC,y=reorder(ID,-log2FC),fill=Col)) + geom_bar() + theme_bw()
heat.in <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[selSpecies,group %>% dplyr::filter(Group != "NA") %>% rownames()]
speciesoriSTAT <- heat.in
mid <- pheatmap::pheatmap(heat.in,scale = "row")
heat.in <- t(scale(t(heat.in),center=T,scale=T))
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
p1 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(vivarium2 == "WT") %>% rownames()],  column_title = "WT", col = col_fun,cluster_rows = F,show_column_names = F)
p2 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(vivarium2 == "KO") %>% rownames()],  column_title = "KO", col = col_fun,cluster_rows = F,show_column_names = F)
p1+p2
######### genus level #########
selSpecies <- read.table("03.signFeatures/Wilcoxon.g.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.WTvsKO < 0.05) %>% dplyr::filter(WT_Median > 0 | KO_Median > 0) %>% rownames()
read.table("03.signFeatures/Wilcoxon.g.xls",sep = "\t",header = T,row.names = 1) %>% dplyr::filter(wilcox.test.p.WTvsKO < 0.05) %>% dplyr::filter(WT_Median > 0 | KO_Median > 0) %>% dplyr::mutate(log2FC=ifelse(KO_Mean > WT_Mean,log2(KO_Mean/WT_Mean),-log2(WT_Mean/KO_Mean))) %>% tibble::rownames_to_column('ID') %>%  dplyr::mutate(Col=ifelse(log2FC > 0,'red','blue')) %>% ggplot(aes(weight=log2FC,y=reorder(ID,-log2FC),fill=Col)) + geom_bar() + theme_bw()
heat.in <- read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep = "\t",header = T,row.names = 1)[selSpecies,group %>% dplyr::filter(Group != "NA") %>% rownames()]
speciesoriSTAT <- heat.in
mid <- pheatmap::pheatmap(heat.in,scale = "row")
heat.in <- t(scale(t(heat.in),center=T,scale=T))
col_fun = colorRamp2(c(-3, 0, 3), c("#6975ab", "white", "firebrick3"))
p1 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(vivarium2 == "WT") %>% rownames()],  column_title = "WT", col = col_fun,cluster_rows = F,show_column_names = F)
p2 <- Heatmap(heat.in[mid$tree_row$labels[mid$tree_row$order],group %>% dplyr::filter(vivarium2 == "KO") %>% rownames()],  column_title = "KO", col = col_fun,cluster_rows = F,show_column_names = F)
p1+p2