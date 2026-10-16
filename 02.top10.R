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
############################## top 10 tax ######################
total_color <- c(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5'),pal_aaas()(10),pal_jco()(10),rev(pal_igv()(10)))
for (rank in c('p','g','s')) {
  for (name in analysis.list) {
    mid.name <- NULL
    plotData <- list()
    orderData <- list()
    taxOrder <- list()
    plotRe <- list()
    for (mid.lv in analysis.list.info$levels[[name]]) {
      dat1 <- read.delim(paste0("00.RawData/",rank,'_tax_rel.sort.screen.txt'), row.names = 1, sep = '\t', head = TRUE, check.names = FALSE)[,group %>% dplyr::filter(get(name) == mid.lv) %>% rownames()]
      dat1 <- dat1[,colSums(dat1) > 0]
      dat1 <- dat1[names(sort(apply(dat1,1,mean),decreasing=T)),] %>% head(n=15)
      dat1 <- dat1[,colSums(dat1) > 0]
      dis_bray <- vegan::vegdist(t(dat1), method = 'bray')
      tree1 <- hclust(dis_bray, method = 'average')
      dat1.plot <- dat1 %>% tibble::rownames_to_column(var="Tax") %>% reshape2::melt(by="Tax") 
      dat1.plot$Tax <- factor(dat1.plot$Tax,levels = rownames(dat1) %>% rev)
      mid.name <- c(mid.name,rownames(dat1))
      plotData[[mid.lv]] <- dat1.plot
      orderData[[mid.lv]] <- colnames(dat1)[tree1$order]
      taxOrder[[mid.lv]] <- rownames(dat1)
    }
    mid.name <- (mid.name %>% table %>% as.data.frame() %>% dplyr::arrange(-Freq))[,1] %>% as.vector()
    phylum_color <- total_color[1:length(mid.name)]
    names(phylum_color) <- mid.name
    
    ggplot(plotData[[mid.lv]],aes(variable,weight=value,fill=Tax)) + geom_bar() + xlim(orderData[[mid.lv]]) + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust = 1),axis.ticks.x = element_blank()) + scale_fill_manual(values = phylum_color[taxOrder[[mid.lv]]]) + xlab(mid.lv) + ylab("Relative Abundance") + guides(fill = guide_legend( ncol = 1, byrow = F))
    
  }
}