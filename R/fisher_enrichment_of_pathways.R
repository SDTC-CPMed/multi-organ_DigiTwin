library(readxl)
library(dplyr)

rm(list = ls())


enrichment <- function(AIDs, group){
  AIDs$X.log.p.value. <- gsub('[,]', '.', AIDs$X.log.p.value.)
  AIDs$X.log.p.value. <- gsub('−', '-', AIDs$X.log.p.value.)
  AIDs$X.log.p.value. <- gsub('×10', 'e', AIDs$X.log.p.value.)
  AIDs$X.log.p.value. <- gsub('\\^', '', AIDs$X.log.p.value.)
  AIDs$X.log.p.value. <- as.numeric(AIDs$X.log.p.value.)
  AIDs <- AIDs[which(AIDs$X.log.p.value. > -log10(0.05)),]
  
  CIA_subs <- unique(sort(CIA$subclusters))
  CIA_mains <- unique(sort(CIA$clustersGlobal))
  
  output_enrichements_pvals <- matrix(NA, ncol = length(CIA_subs) + length(CIA_mains), nrow = 1)
  colnames(output_enrichements_pvals) <- c(CIA_mains, CIA_subs)
  rownames(output_enrichements_pvals) <- group
  
  output_enrichements_odds <- matrix(NA, ncol = length(CIA_subs) + length(CIA_mains), nrow = 1)
  colnames(output_enrichements_odds) <- c(CIA_mains, CIA_subs)
  rownames(output_enrichements_odds) <- group
  
  # i <- 1
  for (i in 1:length(colnames(output_enrichements_pvals))){
    CIA_sub <- colnames(output_enrichements_pvals)[i]
    
    AID.x <- AIDs
    if (length(grep('\\.', CIA_sub)) == 0){
      CIA.x <- CIA[CIA$clustersGlobal == CIA_sub,]
    } else {
      CIA.x <- CIA[CIA$subclusters == CIA_sub,]
    }
    
    input_fisher <- matrix(NA, nrow = 2, ncol = 2)
    colnames(input_fisher) <- c('in CIA.x', 'not in CIA.x')
    rownames(input_fisher) <- c('in AID.x', 'not in AID.x')
    input_fisher[1,1] <- length(which(CIA.x$IngenuityCanonicalPathways %in% AID.x$Ingenuity.Canonical.Pathways))
    input_fisher[1,2] <- length(AID.x$Ingenuity.Canonical.Pathways) - input_fisher[1,1]
    input_fisher[2,1] <- length(CIA.x$IngenuityCanonicalPathways) - input_fisher[1,1]
    
    
    bg_CIA <- CIA
    bg_AID <- AIDs
    
    bg_paths <- unique(sort(c(bg_CIA$IngenuityCanonicalPathways, bg_AID$Ingenuity.Canonical.Pathways)))
    bg_paths_notin <- bg_paths[which(bg_paths %in% CIA.x$IngenuityCanonicalPathways == F)]
    bg_paths_notin <- bg_paths_notin[which(bg_paths_notin %in% AID.x$Ingenuity.Canonical.Pathways == F)]
    
    
    input_fisher[2,2] <- length(bg_paths_notin)
    
    output_fisher <- fisher.test(input_fisher, alternative = 'greater')
    
    output_enrichements_odds[1,i] <- as.data.frame(output_fisher$estimate)[1,1]
    output_enrichements_pvals[1,i] <- output_fisher$p.value
  }
  output_enrichements_fdr <- matrix(p.adjust(output_enrichements_pvals, method = 'fdr'), ncol = length(colnames(output_enrichements_pvals)), nrow = length(rownames(output_enrichements_pvals)))
  colnames(output_enrichements_fdr) <- colnames(output_enrichements_pvals)
  rownames(output_enrichements_fdr) <- rownames(output_enrichements_pvals)
  return(list(fdr = output_enrichements_fdr, odds = output_enrichements_odds, pvals = output_enrichements_pvals))
  
}


# load the data
CIA <- read.csv('../data/Connective_pathway_analysis/TreeStructure_nodes2_AID_noblood.txt', sep = '\t') # all IMIDs
all_CIA <- read.csv('../data/Connective_pathway_analysis/TreeStructure_nodes2_CLUSTER2_AID_noblood.txt', sep = '\t') # UC
CIA <- full_join(CIA, all_CIA)
remove(all_CIA)

# sort out relevant columns
CIA$subclusters <- paste('IMIDs_', CIA$clustersGlobal, '.', CIA$subclusters, sep = '') 
CIA$clustersGlobal <- paste('IMIDs_', CIA$clustersGlobal, sep = '') 
CIA <- CIA[,which(colnames(CIA) %in% c('IngenuityCanonicalPathways', 'clustersGlobal', 'subclusters'))]


AIDs_responders <- read.csv('../data/GSE92415/Gse92415UC_untreate responder vs control_pathway.txt', sep = '\t', skip = 1)
  
res_responders = enrichment(AIDs_responders, 'responders')

AIDs_nonresponders <- read.csv('../data/GSE92415/Gse92415UC_untreated non responder vs control_pathway.txt', sep = '\t', skip = 1)

res_nonresponders = enrichment(AIDs_nonresponders, 'nonresponders')

res = data.frame(as.double(c(res_responders$fdr, res_nonresponders$fdr)))
res[,2] = c(colnames(res_responders$fdr), colnames(res_nonresponders$fdr))
res[,3] = c(rep('responders', 25), rep('nonresponders', 25))
res[,4] = as.factor(res[,1] < 0.05)
colnames(res) = c('fdr', 'program', 'group', 'significant')

#ggplot(res,aes(x=state, y=program, size = fdr))+
#  geom_point()

library(ggplot2)
temp_plot = ggplot(res,aes(x=group, y=program, color=significant, size = -log10(fdr))) +
  geom_point() +
  scale_size(range = c(0, 12),name="-log10(BH fdr P value)",guide = "legend",limits = c(0,60))+
  theme(axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size=13), #change legend title font size
      legend.text = element_text(size=13))
temp_plot

#ggsave(temp_plot, file=paste('/Users/danga10/Documents/CIA/DS_UR_Program_test/tree_structures_code_v2/overlap_diseases_programs/',
#                             'GWAS_overlap_diseases_programs',disease_name,'.pdf',sep=''),
#       width = 8, height = 8)

