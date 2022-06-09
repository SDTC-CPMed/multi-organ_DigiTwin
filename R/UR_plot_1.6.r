library(ggplot2)


URs = read.csv('../data/UR_analysis/UR_IMID_summary_z.csv')

res = data.frame(rep(URs$X0, 32))
res[,2] = as.numeric(as.vector(as.matrix(URs[,2:33])))
res[,3] = as.vector(t(matrix(c(rep(colnames(URs)[2:33], 8)), ncol = 8)))
res[,4] = res[,2] > 0
colnames(res) = c('UR', 'z_score', 'dataset', 'positive')
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

temp_plot = ggplot(res,aes(x=UR, y=dataset, color=positive, size = abs(z_score))) +
  geom_point() +
  geom_hline(yintercept = 23.5) +
  scale_size(range = c(0, 6),name="z_score",guide = "legend",limits = c(0,60))+
  theme(axis.text.x = element_text(angle = 45, size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20),
        legend.title = element_text(size=13), #change legend title font size
        legend.text = element_text(size=13),
        #legend.position = "bottom",
        #legend.justification = "right",
        legend.box = "vertical")+
  #guides(size = guide_legend(nrow = 1), color = guide_legend(nrow =1)) +
  #guides(size = FALSE, color = FALSE) +
  colScale
temp_plot

URs = read.csv('../data/UR_analysis/UR_IMID_summary_logFC.csv')

res = data.frame(rep(URs$X0, 32))
res[,2] = as.numeric(as.vector(as.matrix(URs[,2:33])))
res[,3] = as.vector(t(matrix(c(rep(colnames(URs)[2:33], 8)), ncol = 8)))
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
        #legend.position = "bottom",
        #legend.justification = "right",
        legend.box = "vertical")+
  #guides(size = guide_legend(nrow = 1), color = guide_legend(nrow =1)) +
  #guides(size = FALSE, color = FALSE) +
  colScale
temp_plot






overlap = read.csv('../data/IMIDs_pathway_overlap_with_SPs_reshaped_for_dot_plotALL.txt')

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
  #ggtitle('IMIDs pathway overlap') +
  theme(axis.text.x = element_text(angle = 45, size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20),
        legend.title = element_text(size=13), #change legend title font size
        legend.text = element_text(size=13),
        legend.position = "bottom",
        legend.justification = "right",
        legend.box = "horizontal")+
  guides(size = guide_legend(nrow = 1), color = guide_legend(nrow =1)) +
  colScale
temp_plot





###################################



GWASP2 = read.table('../data/Pval_summary_AID_cluster2_oneSided.csv',header=T,sep=';',dec = ",")
#saveName ='P2'
hv = seq(1:(ncol(GWASP2)-1))
hv = paste('2.',as.character(hv),sep='')
names(GWASP2) = c('dis_id',hv)
GWASP1 = read.table('../data/Pval_summary_AID_cluster1_oneSided.csv',header=T,sep=';',dec = ",")
#saveName ='P1'
hv = seq(1:(ncol(GWASP1)-1))
hv = paste('1.',as.character(hv),sep='')
names(GWASP1) = c('dis_id',hv)
sum(rownames(GWASP1)!=rownames(GWASP2))
GWAS = cbind(GWASP1,GWASP2[,2:ncol(GWASP2)])


GWAS$active.dis = rep('inflamed',length(GWAS$dis_id))
GWAS$active.dis[grep('inactive',GWAS$dis_id)]='non-inflamed'
GWAS$active.dis[GWAS$dis_id=='GSE66413_at_risk_T1D_pancreatic_lymph_nodes']='non-inflamed'
sum(GWAS$active.dis=='non-inflamed')
GWAS$dis.category = GWAS$dis_id
GWAS$dis.category[grep('_UC_',GWAS$dis.category)]='UC'
GWAS$dis.category[grep('_CD_',GWAS$dis.category)]='CD'
GWAS$dis.category[grep('_lupus_',GWAS$dis.category)]='lupus'
GWAS$dis.category[grep('_DLE_',GWAS$dis.category)]='lupus'
GWAS$dis.category[grep('_SLE_',GWAS$dis.category)]='lupus'
GWAS$dis.category[grep('_LN_',GWAS$dis.category)]='lupus'
GWAS$dis.category[grep('_SCLE_',GWAS$dis.category)]='lupus'
GWAS$dis.category[grep('_cSLE_',GWAS$dis.category)]='lupus'
GWAS$dis.category[grep('_AD_',GWAS$dis.category)]='AD'
GWAS$dis.category[grep('_SS_',GWAS$dis.category)]='SS'
GWAS$dis.category[grep('_JM_',GWAS$dis.category)]='JM'
GWAS$dis.category[grep('_PSO_',GWAS$dis.category)]='PSO'
GWAS$dis.category[grep('_RA_',GWAS$dis.category)]='RA'
GWAS$dis.category[grep('_SSc_',GWAS$dis.category)]='SSc'
GWAS$dis.category[grep('_at_risk_T1D_',GWAS$dis.category)]='at_risk_T1D'

clustpos = colnames(GWAS)[2:24]
GWAS$dis.category[GWAS$active.dis=='inflamed']=paste(GWAS$dis.category[GWAS$active.dis=='inflamed'],'_inflamed', sep = '')
GWAS$dis.category[GWAS$active.dis=='non-inflamed']=paste(GWAS$dis.category[GWAS$active.dis=='non-inflamed'],'_non-inflamed', sep = '')
udis = unique(GWAS$dis.category)
P.combined = matrix(0,length(udis),length(clustpos))
rownames(P.combined) = udis
colnames(P.combined) = colnames(GWAS[,clustpos])

idx.p.combined=1
for (z in 1:length(udis)){
  idx.dis = GWAS$dis.category==udis[z]
  for (zz in clustpos){
    Chi2 = -2*sum(log(GWAS[idx.dis,zz]))
    P.combined[z,idx.p.combined] = 1-pchisq(Chi2, 2*length(GWAS[idx.dis,zz]))
    idx.p.combined=idx.p.combined+1
    rm('Chi2')
  }
  idx.p.combined=1
}


d = P.combined[,'1.6']

d[grepl('_inflamed',names(d))]=d[grepl('_inflamed',names(d))]+100
idx = sort(d,decreasing=T,index.return = T)
levels.dis = names(idx$x)

d=P.combined
d[d>0.05]=1
d[d==0]=10^-20
d= -log10(d)
idx = sort(colSums(d),decreasing=T,index.return = T)
levels.programs = colnames(P.combined[,idx$ix])


P.combined.reshaped = data.frame(dis = rep(rownames(P.combined),times = ncol(P.combined)), 
                                 sub.program = rep(colnames(P.combined),each=nrow(P.combined)),
                                 P = as.vector(as.matrix(P.combined)))
P.combined.reshaped$qval=p.adjust(P.combined.reshaped$P, method = "BH", n = length(P.combined.reshaped$P))

#P.combined.reshaped$minus.log10.P = -log10(P.combined.reshaped$P)
P.combined.reshaped$minus.log10.qval =P.combined.reshaped$qval
P.combined.reshaped$minus.log10.qval[P.combined.reshaped$minus.log10.qval==0]=10^-20
P.combined.reshaped$minus.log10.qval = -log10(P.combined.reshaped$minus.log10.qval)

P.combined.reshaped$significant  = P.combined.reshaped$minus.log10.qval>-log10(0.05)
P.combined.reshaped$dis = factor(P.combined.reshaped$dis,levels = levels.dis)

for(i in 1:345){
  P.combined.reshaped$sub.program[i] = paste('IMID_SP',P.combined.reshaped$sub.program[i], sep = '')
}
#P.combined.reshaped$sub.program = factor(P.combined.reshaped$sub.program,levels = levels.programs)

P_order = c(6, 16, 10, 2, 15, 5, 17, 1, 9, 19, 4, 8, 12, 3, 7, 18, 23, 21, 11, 13, 14, 20, 22)
disease_order = c(9, 4, 11, 7, 15, 2, 3, 8, 12, 5, 14, 13, 1, 10, 6)

order = rep(0, 345)
k = 1
for (i in 1:length(P_order)){
  for (j in 1:length(disease_order)){
    order[k] = (P_order[i]-1) * length(disease_order) + disease_order[j]
    k = k+1
  }
}
P.combined.reshaped = P.combined.reshaped[order,]



P.combined.reshaped$sub.program = factor(P.combined.reshaped$sub.program,
                                         levels=unique(P.combined.reshaped$sub.program))
P.combined.reshaped$dis = factor(P.combined.reshaped$dis, levels=unique(P.combined.reshaped$dis))

myColors <- c('red', '#00A9FF')
names(myColors) <- levels(P.combined.reshaped$significant)
colScale <- scale_colour_manual(name = "significant",values = myColors)

temp_plot = ggplot(P.combined.reshaped, aes(x=sub.program, y=dis, color=significant, size = minus.log10.qval)) +
  geom_point() +
  geom_hline(yintercept = 9.5) +
  scale_size(range = c(0, 12),name="-log(P value)",guide = "legend",limits = c(0,100))+
  #ggtitle('GWAS enrichment') +
  theme(axis.text.x = element_text(angle = 45, size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20),
        legend.title = element_text(size=13), #change legend title font size
        legend.text = element_text(size=13),
        legend.position = "bottom",
        legend.justification = "right",
        legend.box = "horizontal")+
  
  theme(text=element_text(size=16, 
                          #       family="Comic Sans MS"))
                          #       family="CM Roman"))
                          #       family="TT Times New Roman"))
                          #       family="Sans"))
                          family="Serif"))


  guides(size = guide_legend(nrow = 1), color = guide_legend(nrow =1)) +
  colScale
temp_plot
