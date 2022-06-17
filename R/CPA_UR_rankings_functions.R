runCPA = function(xx,pathinfo,SAVENAME,clusterscut.subprograms){
  # parts of the code written based on Fredrik Barrenas code
  require('pheatmap')
  
  if (!exists('clusterscut.subprograms')){
    clusterscut.subprograms=1.78
  }
  # remove everything that sits on the diagonal:
  for(ii in 1:nrow(xx)){xx[ii,ii] = NA}
  rownames(xx)= colnames(xx)
  
  # add rownames to the pathinfo
  rownames(pathinfo) = pathinfo$key
  
  
  # cut the tree to find main programs:
  hc.x = hclust(as.dist(1-xx),method='ward.D2') # xx contains Jaccard Index - symmetric matrix same ordering of pathways in columns and in rows
  clustersGlobal = cutree(hc.x,k=2)
  clustersGlobal = as.data.frame(clustersGlobal)
  clustersGlobal$key = rownames(clustersGlobal)
  pathinfo = merge(pathinfo,clustersGlobal, by='key')
 
  ## cut program 1 to find sub-programs:
  for( main.program in unique(pathinfo$clustersGlobal)){
    v.selected.label = clustersGlobal$key[clustersGlobal$clustersGlobal==main.program]
    xx.selected = xx[v.selected.label,v.selected.label]
    pathinfo.selected = pathinfo[v.selected.label,]
    hc.x = hclust(as.dist(1-xx.selected),method='ward.D2')
    subclusters = cutree(hc.x, h=clusterscut.subprograms)
    # table(subclusters)
    
    #subclusters = cutree(hc.x,k=10)
    subclusters = as.data.frame(subclusters)
    subclusters$key = rownames(subclusters)
    subclusters$subclusters = paste(main.program,subclusters$subclusters,sep='.')
    if(!exists('subclusters.collected')){
      subclusters.collected = subclusters
    }else{
      subclusters.collected = rbind(subclusters.collected,subclusters)
    }
    
  }
 
  pathinfo = merge(pathinfo,subclusters.collected, by='key')
  return(pathinfo)
  
}

count.same.and.opposing.activations = function(pathinfo){
  
  pathinfo$opposit.direction = rep(0,length(pathinfo$key))
  idx = (pathinfo$inactiveTissue=='inhibited' & pathinfo$activeTissue == 'activated') | 
    (pathinfo$inactiveTissue=='activated' & pathinfo$activeTissue == 'inhibited') |
    (pathinfo$inactiveTissue=='not_sig'   & pathinfo$activeTissue != 'not_sig')   |
    (pathinfo$inactiveTissue!='not_sig'   & pathinfo$activeTissue == 'not_sig')
  pathinfo$opposit.direction[idx] = 1
  
  rm(list='idx')
  pathinfo$same.direction = rep(0,length(pathinfo$key))
  idx = (pathinfo$inactiveTissue=='activated' & pathinfo$activeTissue == 'activated') | (pathinfo$inactiveTissue=='inhibited' & pathinfo$activeTissue == 'inhibited') #|
  pathinfo$same.direction[idx] = 1
  
  ## here I check the percentage of pathways showing opposit activation profile:
  u.clusters = data.frame(cluster.id = unique(pathinfo$clustersGlobal), ratio = rep(NA,length(unique(pathinfo$clustersGlobal))))
  for (z in u.clusters$cluster.id){
    u.clusters$ratio[z] = mean(pathinfo$opposit.direction[pathinfo$clustersGlobal==z])
  }
  idx = sort(u.clusters$ratio,decreasing=T,index.return =T)
  u.clusters.main = u.clusters[idx$ix,]

  u.clusters = data.frame(cluster.id = unique(pathinfo$subclusters), ratio = rep(NA,length(unique(pathinfo$subclusters))))
  for (z in 1:length(u.clusters$cluster.id)){
    u.clusters$ratio[z] = mean(pathinfo$opposit.direction[pathinfo$subclusters==u.clusters$cluster.id[z]])
  }
  u.clusters = rbind(u.clusters.main,u.clusters)
  names(u.clusters) = c('Cluster','Percentage.opposit.direction')
  u.clusters$Percentage.opposit.direction = round(u.clusters$Percentage.opposit.direction*100)
  
  
  
  ## here I check the percentage of pathways showing same activation profile:
  u.clusters2 = data.frame(cluster.id = unique(pathinfo$clustersGlobal), ratio = rep(NA,length(unique(pathinfo$clustersGlobal))))
  for (z in u.clusters2$cluster.id){
    u.clusters2$ratio[z] = mean(pathinfo$same.direction[pathinfo$clustersGlobal==z])
  }
  u.clusters2.main=u.clusters2
  u.clusters2 = data.frame(cluster.id = unique(pathinfo$subclusters), ratio = rep(NA,length(unique(pathinfo$subclusters))))
  for (z in 1:length(u.clusters2$cluster.id)){
    u.clusters2$ratio[z] = mean(pathinfo$same.direction[pathinfo$subclusters==u.clusters2$cluster.id[z]])
  }
  u.clusters2 = rbind(u.clusters2,u.clusters2.main)
  names(u.clusters2) = c('Cluster','Percentage.same.direction')
  u.clusters2$Percentage.same.direction = round(u.clusters2$Percentage.same.direction*100)
  
  u.clusters = merge(u.clusters,u.clusters2,by = 'Cluster')
  u.clusters = u.clusters[order(u.clusters$Cluster),]
  
  return(u.clusters)
}

FindOverlapingPaths = function(Dis.Pval,pathinfo){
  require('fdrtool')
  udis = unique(Dis.Pval$disease)
  uprog = c(unique(pathinfo$subclusters),c('P1','P2'))
  
  for (u in udis){
    dhv = Dis.Pval[Dis.Pval$disease==u,]
    dhv = dhv[dhv$pvalue<0.05,]
    # dhv = Dis.Pval$IngenuityCanonicalPathways[Dis.Pval$disease==u]
    dhv = dhv$IngenuityCanonicalPathways[dhv$IngenuityCanonicalPathways %in% pathinfo$IngenuityCanonicalPathways]
    for (p in uprog){
      if(grepl('P',p)){
        phv = pathinfo$IngenuityCanonicalPathways[pathinfo$clustersGlobal==gsub('P','',p)]
      }else{
        phv = pathinfo$IngenuityCanonicalPathways[pathinfo$subclusters==p]
      }
      
      a = sum(dhv %in% phv)
      b = sum(!dhv %in% phv)
      c = sum(!phv %in% dhv)
      d = sum(!pathinfo$IngenuityCanonicalPathways %in% c(phv,dhv))
      hv = fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2),alternative='greater')
      hv = data.frame(program = p,udis_condition = u, p_val = hv$p.value, a = a, b = b, c = c, d = d)
      if (!exists('Dis.Pval.Overlap')){
        Dis.Pval.Overlap = hv
      }else{
        Dis.Pval.Overlap = rbind(Dis.Pval.Overlap,hv)
      }
      rm(list=c('hv','a','b','c','d','phv'))
    }
    rm('dhv')
  }
  Dis.Pval.Overlap$bhfdr = p.adjust(p = Dis.Pval.Overlap$p_val,method = 'BH')
  # hv=fdrtool(x=Dis.Pval.Overlap$p_val,statistic='pvalue',plot=F)
  # Dis.Pval.Overlap$bhfdr = hv$lfdr
  return(Dis.Pval.Overlap)
}