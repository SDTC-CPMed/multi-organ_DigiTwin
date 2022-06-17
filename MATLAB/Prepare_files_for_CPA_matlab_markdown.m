clear all
close all
clc

%% Define main path
InputOutputFiles = '../data/CPA_InputFiles/';

%% CIA scRNAseq:

% define variables:
savename = 'scRNAseq';
PATHH1 = sprintf('%sPathwayenrichment_results/scRNAseq_CIA/',InputOutputFiles);
FN = readtable(sprintf('%s/PathFilesDescription_CIA.csv',InputOutputFiles),'delimiter',',');

% read in all pathways, filter significant ones, find all involved genes
% and summarize pathway activation direction:
AllEnrichedPaths = ReadInIPAPathwayEnrichments(FN,PATHH1);

% calculate jaccard index:
Jaccard = CalculateJaccardIndex(AllEnrichedPaths);

% save to files:
writetable(Jaccard,sprintf('%sJaccardIndex_%s.csv',InputOutputFiles,savename))
writetable(AllEnrichedPaths,sprintf('%sPathInfo_%s.csv',InputOutputFiles,savename))
writetable(FN,sprintf('%sDatasetInfo_%s.csv',InputOutputFiles,savename),'delimiter','\t')


%% IMIDs:

% define variables:
savename = 'IMID';
PATHH1 = sprintf('%sPathwayenrichment_results/IMIDs/',InputOutputFiles);
FN = readtable(sprintf('%s/PathFilesDescription_IMIDs.csv',InputOutputFiles));

% read in all pathways, filter significant ones, find all involved genes
% and summarize pathway activation direction:
[AllEnrichedPaths,filesForOverlap] = ReadInIPAPathwayEnrichments(FN,PATHH1);

% calculate jaccard index:
Jaccard = CalculateJaccardIndex(AllEnrichedPaths);

% save to files:
writetable(Jaccard,sprintf('%sJaccardIndex_%s.csv',InputOutputFiles,savename))
writetable(AllEnrichedPaths,sprintf('%sPathInfo_%s.csv',InputOutputFiles,savename))
writetable(FN,sprintf('%sDatasetInfo_%s.csv',InputOutputFiles,savename),'delimiter','\t')
writetable(filesForOverlap,sprintf('%sDieasesPathways_%s.csv',InputOutputFiles,savename),'delimiter','\t')

%% RUN R script to get CPA done

