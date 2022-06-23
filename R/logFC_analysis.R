logFC_analysis <- function(Datasets_Inf, Datasets_Noninf, path_DEGs){
  i = 0
  logFC_list = list()
  for dataset in Datasets_Inf.columns:
    data = pd.read_table(path_DEGs + dataset)
  if i in (6,12,15,16):
    if i == 15:
    data = data.rename(columns = {'ORF': 'ENTREZ_GENE_ID'})
  trans = translation[translation['Symbol'].isin(Datasets_Inf.index)][['GeneID','Symbol']]
  data = trans.merge(data, left_on = 'GeneID', right_on ='ENTREZ_GENE_ID')[['Symbol', 'logFC']]
  data = data.rename(columns = {'Symbol': 'Gene.symbol'})
  else:
    if i in (7,20):
    data = data.rename(columns = {'ORF': 'Gene.symbol'})
  if i == 8:
    data = data.rename(columns = {'Gene.Symbol': 'Gene.symbol'})
  if i == 9:
    data = data.rename(columns = {'GENE_SYMBOL': 'Gene.symbol'})
  if i == 13:
    data = data.rename(columns = {'ID': 'Gene.symbol'})
  data = data[data['adj.P.Val'] < 0.05]
  data = data[data['Gene.symbol'].isin(Datasets_Inf.index)]\
  .drop_duplicates('Gene.symbol')[['Gene.symbol', 'logFC']]
  
  data = data.rename(columns = {'logFC': dataset})
  data.index = data['Gene.symbol']
  data.pop('Gene.symbol')    
  logFC_list.append(data)
  i = i+1
  
  logFC_Inf = pd.concat(logFC_list, axis = 1)
  
  i = 0
  logFC_list = list()
  for dataset in Datasets_Noninf.columns:
    
    data = pd.read_table(path_DEGs + dataset)
  if i == 3:
    data = data.rename(columns = {'ORF': 'Gene.symbol'})
  if i == 8:
    data = data.rename(columns = {'GENE_NAME': 'Gene.symbol'})
  data = data[data['adj.P.Val'] < 0.05]
  data = data[data['Gene.symbol'].isin(Datasets_Inf.index)].drop_duplicates('Gene.symbol')[['Gene.symbol', 'logFC']]
  data = data.rename(columns = {'logFC': dataset})
  data.index = data['Gene.symbol']
  data.pop('Gene.symbol')
  logFC_list.append(data)
  i = i+1
  
  logFC_Noninf = pd.concat(logFC_list, axis = 1)
  
  return(list(logFC_Inf, logFC_Noninf))
}