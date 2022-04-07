library(ggplot2)


URs = read.csv('../data/UR_analysis/UR_IMID_summary_z.csv')

res = data.frame(rep(URs$X0, 32))
res[,2] = as.numeric(as.vector(as.matrix(URs[,2:33])))
res[,3] = as.vector(t(matrix(c(rep(colnames(URs)[2:33], 8)), ncol = 8)))
res[,4] = res[,2] > 0
colnames(res) = c('UR', 'z_score', 'dataset', 'positive')
res$dataset = factor(res$dataset, levels=unique(res$dataset))


temp_plot = ggplot(res,aes(x=UR, y=dataset, color=positive, size = abs(z_score))) +
  geom_point() +
  geom_hline(yintercept = 23.5) +
  scale_size(range = c(0, 12),name="z_score",guide = "legend",limits = c(0,60))+
  theme(axis.text.x = element_text(angle = 45))
temp_plot

URs = read.csv('../data/UR_analysis/UR_IMID_summary_logFC.csv')

res = data.frame(rep(URs$X0, 32))
res[,2] = as.numeric(as.vector(as.matrix(URs[,2:33])))
res[,3] = as.vector(t(matrix(c(rep(colnames(URs)[2:33], 8)), ncol = 8)))
res[,4] = res[,2] > 0
colnames(res) = c('UR', 'logFC', 'dataset', 'positive')
res$dataset = factor(res$dataset, levels=unique(res$dataset))


temp_plot = ggplot(res,aes(x=UR, y=dataset, color=positive, size = abs(logFC))) +
  geom_point() +
  geom_hline(yintercept = 23.5) +
  scale_size(range = c(0, 12),name="logFC",guide = "legend",limits = c(0,60))+
  theme(axis.text.x = element_text(angle = 45))
temp_plot


