plot_overlap <- function(overlap){
  overlap[,3] = -log10(overlap[,4])
  overlap[,4] = overlap[,4] < 0.05
  colnames(overlap) = c('program', 'disease', 'bhfdr', 'significant')
  
  for(i in 1:375){
    if(!overlap$program[i] %in% c('P1', 'P2')){
      overlap$program[i] = paste('IMID_SP',overlap$program[i], sep = '')
    }
  }
  
  P_order = c(1, 8, 6, 7, 9, 3, 5, 11, 21, 10, 4, 19, 22, 23, 12, 14, 15, 16, 20, 25, 13, 17, 18, 24, 2)
  disease_order = c(11, 15, 1, 12, 9, 3, 7, 10, 5, 8, 4, 2, 14, 13, 6)
  
  order = rep(0, 375)
  k = 1
  for (i in 1:length(disease_order)){
    for (j in 1:length(P_order)){
      order[k] = (disease_order[i]-1) * length(P_order) + P_order[j]
      k = k+1
    }
  }
  overlap = overlap[order,]
  overlap$program = factor(overlap$program, levels=unique(overlap$program))
  overlap$disease = factor(overlap$disease, levels=unique(overlap$disease))
  
  
  myColors <- c('red', '#00A9FF')
  names(myColors) <- levels(overlap$significant)
  colScale <- scale_colour_manual(name = "significant",values = myColors)
  
  temp_plot = ggplot(overlap, aes(x=program, y=disease, color=significant, size = bhfdr)) +
    geom_point() +
    geom_hline(yintercept = 9.5) +
    scale_size(range = c(0, 12),name="-log(P value)",guide = "legend",limits = c(0,100))+
    theme(axis.text.x = element_text(angle = 45, size = 13,vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 20),
          legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=9),
          legend.position = "bottom",
          legend.justification = "right",
          legend.box = "horizontal")+
    guides(size = guide_legend(nrow = 1), color = guide_legend(nrow =1)) +
    colScale
  return(temp_plot)
}