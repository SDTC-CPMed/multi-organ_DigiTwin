## Meta analysis of 11 IMIDs for Connective pathway analysis and UR prioritization

### Differential expression analysis

### Connective pathway analysis

Similarly as described for the [CIA analysis - Connective Pathway
Analysis](./CIA-analyses.md), we have performed connective pathway
analysis for the IMID datasets.

``` r
# load pre-computed Jaccard index matrix and pathway information plus define name which is to be used to save results:
xx = as.matrix(read.csv(paste(MainPath, 'JaccardIndex_IMID.csv',sep='')))
pathinfo = read.table(paste(MainPath,'PathInfo_IMID.csv',sep=''),sep=',',header=1) 
SAVENAME = 'IMID'

pathinfo = runCPA(xx,pathinfo,SAVENAME,1.78)
# save results to file:
write.table(file=paste(MainPath,'CPA_',SAVENAME,'_all_pathways.txt',sep=''),pathinfo,sep="\t", col.names=T, row.names=F, quote=F)

# count the ratios of pathways showing same and opposing activation pattern in inflamed vs non-inflamed organ:
ratios = count.same.and.opposing.activations(pathinfo)
```

## Overlap of the CPA of IMIDs vs individual IMIDs

In order to get a better overview on which programs and sub-programs of
CPA are enriched in indyvidual diseases, we calulated a Fisher Exact
Test. To test if the programs (IMID_P1 and IMID_P2) derived from all
analyzed IMIDs overlapped with programs from individual IMIDs in
inflamed and non-inflamed organ sites, separately, we performed Fisher’s
exact tests (right tailed), using all pathways in connective pathway
analysis as a background, followed by correction for multiple testing
using the Benjamini-Hochberg procedure. These analyses were repeated for
IMID_subprograms (IMID_SPs) and subprograms from individual IMIDs.
Disease pathways were defined as all pathways significantly enriched in
a particular disease and inflammation state (IPA, p \< 0.05). In case
two or more datasets were representative of the same disease and
condition (for example UC inflamed organs), pathway enrichment p values
were combined with Fisher’s method. Enriched pathways were considered
those whose combined p \< 0.05

``` r
# prepare files to check overlap between programs and enriched pathways per IMID:
Dis.Pval = read.table(paste(MainPath,'DieasesPathways_IMID.csv',sep=''),sep='\t',quote = "",header=T)
OverlapingPaths = FindOverlapingPaths(Dis.Pval,pathinfo)
```
