plot_logFC <- function(URs_logFC){
  res = data.frame(rep(URs_logFC$X0, 32))
  res[,2] = as.numeric(as.vector(as.matrix(URs_logFC[,2:33])))
  res[,3] = as.vector(t(matrix(c(rep(colnames(URs_logFC)[2:33], 8)), ncol = 8)))
  res[,4] = res[,2] > 0
  colnames(res) = c('UR', 'logFC', 'dataset', 'positive')
  res$dataset = factor(res$dataset, levels=unique(res$dataset))
  
  myColors <- c('red', '#00A9FF')
  names(myColors) <- levels(res$positive)
  colScale <- scale_colour_manual(name = "positive",values = myColors)
  
  new_names = c('GSE16161_AD', 'GSE32924_AD_infl.', 'GSE16879_CD-colon', 'GSE179285_CD-colon', 'GSE16879_CD-ileum',
                'GSE179285_CD-ileum', 'GSE81071_DLE', 'GSE148810_JM_infl.', 'GSE32591_LN', 'GSE181318_Pso', 'GSE1919_RA',
                'GSE55235_RA', 'GSE81071_SCLE', 'GSE176510_SS', 'GSE40568_SS', 'GSE81292_SSc', 'GSE95065_SSc',
                'GSE11223_UC-desc.colon', 'GSE11223_UC-sigm.colon', 'GSE179285_UC_infl.', 'GSE148810_cSLE', 'GSE112943_lupus',
                'GSE112943_SCLE', 'GSE32924_AD', 'GSE179285_CD', 'GSE75214_CD', 'GSE148810_JM', 'GSE14905_Pso', 
                'GSE75214_UC', 'GSE11223_UC', 'GSE179285_UC', 'GSE66413_T1D')
  
  levels(res$dataset) = new_names
  
  temp_plot = ggplot(res,aes(x=UR, y=dataset, color=positive, size = abs(logFC))) +
    geom_point() +
    geom_hline(yintercept = 23.5) +
    scale_size(range = c(0, 12),name="logFC",guide = "legend",limits = c(0,40))+
    theme(axis.text.x = element_text(angle = 45, size = 13,vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 20),
          legend.title = element_text(size=13), #change legend title font size
          legend.text = element_text(size=13),
          legend.box = "vertical")+
    #theme_bw()+
    colScale
  return(temp_plot)
}