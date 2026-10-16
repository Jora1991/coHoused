
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

############################## diversity ###############################
# plots for alpha diversity
alpha <- read.table("00.rawdata/alpha_diversity_index.xls",sep="\t",header = T,row.names = 1)
alpha <- merge(group %>% dplyr::select(analysis.list) %>% as.data.frame(),alpha,by="row.names")
alpha.re <- matrix(nrow = 4*nrow(alpha),ncol=2+length(analysis.list))
n <- 1
for (i in c('observed_otus','chao1','shannon','simpson')) {
  for (j in 1:nrow(alpha)) {
    alpha.re[n,] <- c(i,alpha[j,i],as.vector(unlist(alpha[j,analysis.list])))
    n <- n+1
  }
}
colnames(alpha.re) <- c("Factor","Abun",analysis.list)
alpha.re <- data.frame(alpha.re)
alpha.re$Abun <- as.numeric(alpha.re$Abun)
for (name in analysis.list){
  box.in <- alpha.re[,c("Abun","Factor", "Group")]
  colnames(box.in) <- c("Abun","Factor", "Group")
  box.in <- box.in[! is.na(box.in$Group), ]
  box.in$Group <- factor(box.in$Group,levels = analysis.list.info$levels[[name]])
  box.in <- box.in %>% filter(Factor %in% c('observed_otus','chao1','shannon','simpson'))
  box.in$Factor <- factor(box.in$Factor,levels = rev(c('observed_otus','chao1','shannon','simpson')) )
  plot.alpha <-  ggplot(box.in, aes(Group, Abun, fill = Group)) +
    labs(x="",y="")+
    #geom_point(position=position_jitter(0.3),aes(fill=Group),cex=2.5)+ #angle=45,vjust=1,hjust=1,"jitter"
    geom_violin() + geom_boxplot(col='black',size=1,width =0.18,outlier.colour = NA) + 
    scale_fill_manual(values = analysis.list.info$colors[[name]])+ theme_classic() +
    theme(axis.text.x=element_text(colour = "black",size=7,angle = 45,hjust = 1), axis.text.y=element_text(colour = "black",size=7),axis.title.y=element_text(size=7),axis.title.x=element_text(size=7),legend.position = "none",axis.text = element_text(size=7)) +
    facet_wrap(vars(Factor), scales = "free",nrow = 1) +
    geom_signif(comparisons = analysis.list.info$compaired[[name]],
                color="black",step_increase = 0.08,
                map_signif_level=function(p) if(p < 0.01){'+'}else{if(p < 0.05){'*'}else{sprintf("%.2f",p)}},
                test = wilcox.test)
}

data <- read.table("00.rawdata/s_tax_rel.sort.screen.txt",sep="\t",header = T,row.names = 1) 
dist <- vegdist(t(data), method="bray")
write.table(as.matrix(dist),file="00.rawdata/species_bray_dm.txt",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")
data <- read.table("00.rawdata/g_tax_rel.sort.screen.txt",sep="\t",header = T,row.names = 1)
dist <- vegdist(t(data), method="bray")
write.table(as.matrix(dist),file="00.rawdata/genus_bray_dm.txt",sep="\t",col.names = NA, row.names = TRUE,quote = FALSE,fileEncoding = "UTF-8")

# plots for beta diversity
for (name in analysis.list){
  for (rank in c('species',"genus")) {
    group.info <- group[,name]
    names(group.info) <- rownames(group)
    group.info <- group.info[!is.na(group.info)]
    #we dist  
    we.dist <- read.table(paste0("00.rawdata/",rank,"_bray_dm.txt"),sep="\t",header = T,row.names = 1)
    dataTB <- as.dist(we.dist[names(group.info),names(group.info)])
    adonis(dataTB~as.vector(group.info),permutations = 999)
    #PCOA
    PCOA <- pcoa(dataTB, correction="none", rn=NULL)
    result <-PCOA$values[,"Relative_eig"]
    pco1 = as.numeric(sprintf("%.3f",result[1]))*100
    pco2 = as.numeric(sprintf("%.3f",result[2]))*100
    pc = as.data.frame(PCOA$vectors)
    pc$names = rownames(pc)
    xlab=paste("PCoA1 (",pco1,"%)",sep="")
    ylab=paste("PCoA2 (",pco2,"%)",sep="")
    out <- pc[,c(1,2)]
    colnames(out) <- c("PC1","PC2")
    pc$Group <- factor(as.vector(group.info),levels = analysis.list.info$levels[[name]])
    ggplot(pc, aes(Axis.1,Axis.2, fill=Group,color=Group,shape=Group)) +
      labs(x=xlab,y=ylab) +
      geom_hline(yintercept=0,linetype=4,color="grey") +
      geom_vline(xintercept=0,linetype=4,color="grey") +
      geom_point(size=4,alpha=0.7) +
      scale_fill_manual(values = analysis.list.info$colors[[name]])+
      scale_color_manual(values = analysis.list.info$colors[[name]])+
      stat_ellipse(show.legend = F,level = 0.95)+
      theme(axis.text.x=element_text(colour = "black",angle=45,vjust=1,hjust=1,size = 7), axis.text.y=element_text(colour = "black",size = 7),panel.background = element_rect(fill="white",color="black",linetype=1,size=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text = element_text(size=7))
    ggsave(paste("01.Diversity/",name,".",rank,".PCoA.pdf",sep=""),plot.pcoa,device = "pdf",width = 7,height = 6)
  }
}
