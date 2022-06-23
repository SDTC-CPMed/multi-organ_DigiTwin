# UR enrichment analysis of CIA
  

```matlab:Code
clear all
close all
clc
warning('off', 'all')
```

# Input requirements
### For this analysis we will need:

   -  DEGs computed by (TODO link!) 
   -  Orthologous translation 
   -  Results from the connective pathway analysis (TODO link!) 
   -  Predicted URs and their downstreamers from Niche Net (TODO link!) 

### We can find the files at following locations:

\hfill \break


```matlab:Code
% DEGs
path_DEGs = '../data/CIA_mouse_DEGs/';
% translation
path_translation = '../data/Connective_pathway_analysis/orthologous_translation_NicheNet_analysis.txt';
% connective pathway analysis
path_cpa = '../data/Connective_pathway_analysis/TreeStructure_nodes2_scRNAseq.txt';
% URs 
path_URs = '../data/Connective_pathway_analysis/all_outgoing_ligand_activity.txt';
```

### We will produce


A table with combined, FDR corrected pvalues and the count of significant p values across all cell types for each program/subprogram and UR 


### We will save the file here:

\hfill \break


```matlab:Code
path_output = '../data/UR_analysis/UR_predictions_CIA_disease_Pvals.xlsx';
```

  
# Prepare the data
  
### Load DEGs and translate the gene names to human orthologous.

\hfill \break


```matlab:Code

Fn = struct2cell(dir(path_DEGs));
Fn(2:end,:)=[];
Fn(~contains(Fn,'.csv'))=[];
Fn = table(Fn','variablenames',{'fn'});
Fn.Cell = cellfun(@(x) x{5}, regexp(Fn.fn,'\_','split'),'UniformOutput',0);
Fn.Organ = cellfun(@(x) x{6}(1:end-4), regexp(Fn.fn,'\_','split'),'UniformOutput',0);
Fn(~ismember(Fn.Organ,{'Joint','Muscle'}),:)=[];

translation = readtable(path_translation);
for z = 1 : length(Fn.fn)
   hv =readtable(sprintf('%s%s',path_DEGs,Fn.fn{z}));
   hv(~strcmp(hv.is_de_fdr_0_05,'True'),:)=[];
   xx= unique(translation(ismember(translation.Mouse_gene_name,hv.Var1),:));
   Fn.DEG{z} = xx.Gene_name;
   clear hv xx
end

% Restructure the data
Fn.cell_dataset = cellfun(@(x) x{5},regexp(Fn.fn,'\_','split'),'UniformOutput',0);
Fn.Properties.VariableNames(strcmp(Fn.Properties.VariableNames,'Organ'))={'Organ_condition'};
Fn.group = repmat({'Joint'},length(Fn.fn),1);
Fn.group(strcmp(Fn.Organ_condition,'Muscle'))={'Muscle'};
Fn.key = cellfun(@(x,y) sprintf('%s***%s',x,y),Fn.cell_dataset,Fn.group,'UniformOutput',0);

Fn = Fn(:,{'key','cell_dataset','group','DEG'});

head(Fn, 5)
```

| |key|cell_dataset|group|DEG|
|:--:|:--:|:--:|:--:|:--:|
|1|'B-cells***Muscle'|'B-cells'|'Muscle'|2919x1 cell|
|2|'Dendritic-cells***M...|'Dendritic-cells'|'Muscle'|2516x1 cell|
|3|'Endothelial-cells**...|'Endothelial-cells'|'Joint'|1737x1 cell|
|4|'Endothelial-cells**...|'Endothelial-cells'|'Muscle'|2460x1 cell|
|5|'Erythrocytes***Musc...|'Erythrocytes'|'Muscle'|3211x1 cell|

### Load CIA_Ps and CIA_SPs from connective pathway analysis

\hfill \break


```matlab:Code

SPhv = readtable(path_cpa);
SPhv.subprograms = cellfun(@(x) sprintf('1.%s',x), strtrim(cellstr(num2str(SPhv.subclusters))),'UniformOutput',0);
SP = table(unique(SPhv.subprograms),'variablenames',{'SP'});
for z = 1 : length(SP.SP)
    SP.AllMolecules{z} = unique(regexp(strjoin(SPhv.AllMolecules(strcmp(SPhv.subprograms,SP.SP{z})),','),'\,','split'));
end

SP.SP(length(SP.SP)+1)={'P1'};
SP.AllMolecules{end} = unique(regexp(strjoin(SPhv.AllMolecules(SPhv.clustersGlobal==1),','),'\,','split'));
SP.SP(length(SP.SP)+1)={'P2'};
SP.AllMolecules{end} = unique(regexp(strjoin(SPhv.AllMolecules(SPhv.clustersGlobal==2),','),'\,','split'));

head(SP, 5)
```

| |SP|AllMolecules|
|:--:|:--:|:--:|
|1|'1.1'|1x393 cell|
|2|'1.10'|1x108 cell|
|3|'1.2'|1x716 cell|
|4|'1.3'|1x247 cell|
|5|'1.4'|1x351 cell|

  
### load URs from Niche net

\hfill \break


```matlab:Code

UR = readtable(path_URs);
UR(UR.pearson<0,:)=[];
UR.senderOrgan = cellfun(@(x) x{2}, regexp(UR.Sender,'\_','split'),'UniformOutput',0);
UR.targetOrgan = cellfun(@(x) x{2}, regexp(UR.Target,'\_','split'),'UniformOutput',0);
UR(~strcmp(UR.senderOrgan,UR.targetOrgan),:)=[];
UR(~ismember(UR.senderOrgan,{'Muscle','Joint'}),:)=[];
UR.targetcell = cellfun(@(x) x{1}, regexp(UR.Target,'\_','split'),'UniformOutput',0);
UR = unique(UR(:,{'test_ligand','targetcell','targetOrgan','target'}));
UR.Properties.VariableNames = {'UR','cell_dataset','condition','DS'};
UR.DS = regexp(UR.DS,'\/','split');
UR.group(strcmp(UR.condition,'Muscle')) = {'Muscle'};
UR.group(strcmp(UR.condition,'Joint')) = {'Joint'};
UR.key = cellfun(@(x,y) sprintf('%s***%s',x,y),UR.cell_dataset,UR.group,'UniformOutput',0);
UR = UR(:,{'UR','key','cell_dataset','group','DS'});
head(UR, 5)
```

| |UR|key|cell_dataset|group|DS|
|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'APOE'|'B-cells***Muscle'|'B-cells'|'Muscle'|1x64 cell|
|2|'APOE'|'Dendritic-cells***M...|'Dendritic-cells'|'Muscle'|1x50 cell|
|3|'APOE'|'Erythrocytes***Musc...|'Erythrocytes'|'Muscle'|1x69 cell|
|4|'APOE'|'Granulocytes***Musc...|'Granulocytes'|'Muscle'|1x24 cell|
|5|'APOE'|'Macrophages***Muscl...|'Macrophages'|'Muscle'|1x18 cell|

  
# Main analysis 
### Fisher test enrichment of UR DS in SPs:

\hfill \break


```matlab:Code

FPval = cell(length(SP.SP),1);
FOR = cell(length(SP.SP),1);

uUR = unique(UR.UR);
for p = 1 : length(SP.SP)
    for cd = 1 : length(Fn.key)               
            spgenes = SP.AllMolecules{p}(ismember(SP.AllMolecules{p},Fn.DEG{cd}));
            URhv = UR(strcmp(UR.key,Fn.key{cd}),:); %!!!!! URs from more datasets - solved .key now correspond to datasets

            for ur = 1 : length(uUR)
                if sum(strcmp(URhv.UR,uUR{ur}))
                    if length(unique(Fn.DEG{cd})) ~= length(Fn.DEG{cd})
                        error = p;
                    end

                    urds = URhv.DS{strcmp(URhv.UR,uUR{ur})}; 

                    a = sum(ismember(spgenes,urds));
                    b = sum(~ismember(spgenes,urds));
                    c = sum(~ismember(urds,spgenes));
                    d = sum(~ismember(Fn.DEG{cd},[spgenes,urds]));

                    [~, FPval{p}(cd,ur)] = fishertest([a b; c d],'tail','right');
                    FOR{p}(cd,ur) = (a*d)/(b*c);
                else
                    FPval{p}(cd,ur) = 1;
                    FOR{p}(cd,ur) = 0;
                end

                clear urds a b c d
            end
            clear spgenes
    end
end

% Summarize results
FishertestR.cellsdes = SP;
FishertestR.rowdes = Fn;
FishertestR.coldes_URs = uUR;
FishertestR.moreURinfo = UR;
FishertestR.Pval = FPval;
FishertestR.OR = FOR;
FishertestR
```


```text:Output
FishertestR = 
      cellsdes: [12x2 table]
        rowdes: [18x4 table]
    coldes_URs: {31x1 cell}
    moreURinfo: [123x5 table]
          Pval: {12x1 cell}
            OR: {12x1 cell}

```

  
### combine pvalues over joint and muscle:

\hfill \break


```matlab:Code

% For each program and UR, use Fisher's method to combine pvals over all cell types 
for p = 1 : length(FishertestR.Pval)
    for ur = 1 : size(FishertestR.Pval{p},2)
        joint_idx = strcmp(FishertestR.rowdes.group,'Joint');
        chi = -2*sum(log(FishertestR.Pval{p}(joint_idx,ur)));
        CombinedPval_joint(p,ur) = my_chi2cdf(chi,2*sum(joint_idx));
        
        muscle_idx = strcmp(FishertestR.rowdes.group,'Muscle');
        chi = -2*sum(log(FishertestR.Pval{p}(muscle_idx,ur)));
        CombinedPval_muscle(p,ur) = my_chi2cdf(chi,2*sum(muscle_idx));
        
        chi = -2*sum(log(FishertestR.Pval{p}(:,ur)));
        CombinedPval_all(p,ur) = my_chi2cdf(chi,2*32);
    end
end

% Create a summary file where we store important results
CombinedOverJointandMuscle.cols_URs = FishertestR.coldes_URs;
CombinedOverJointandMuscle.rows_SPs = FishertestR.cellsdes.SP;
CombinedOverJointandMuscle.CombP_joint = CombinedPval_joint;
CombinedOverJointandMuscle.CombP_muscle = CombinedPval_muscle;
CombinedOverJointandMuscle.CombP_all = CombinedPval_all;
CombinedOverJointandMuscle
```


```text:Output
CombinedOverJointandMuscle = 
        cols_URs: {31x1 cell}
        rows_SPs: {12x1 cell}
     CombP_joint: [12x31 double]
    CombP_muscle: [12x31 double]
       CombP_all: [12x31 double]

```

  
### FDR correction over Programs

\hfill \break


```matlab:Code

% For each program and UR, do FDR correction 
for i = 1:length(CombinedOverJointandMuscle.CombP_joint(:,1))
    FDR_joint(i,:) = mafdr(CombinedOverJointandMuscle.CombP_joint(i,:), 'BHFDR',true);
    FDR_muscle(i,:) = mafdr(CombinedOverJointandMuscle.CombP_muscle(i,:), 'BHFDR',true);
    FDR_all(i,:) = mafdr(CombinedOverJointandMuscle.CombP_all(i,:), 'BHFDR',true);

end

CombinedOverJointandMuscle.FDR_Joint = FDR_joint;
CombinedOverJointandMuscle.FDR_Muscle = FDR_muscle;
CombinedOverJointandMuscle.FDR_All = FDR_all;

CombinedOverJointandMuscle
```


```text:Output
CombinedOverJointandMuscle = 
        cols_URs: {31x1 cell}
        rows_SPs: {12x1 cell}
     CombP_joint: [12x31 double]
    CombP_muscle: [12x31 double]
       CombP_all: [12x31 double]
       FDR_Joint: [12x31 double]
      FDR_Muscle: [12x31 double]
         FDR_All: [12x31 double]

```

  
  
### Count significant pvalues over joints and muscles:

\hfill \break


```matlab:Code

% For each program and UR, count how in how many cell types was UR significant 
for p = 1 : length(FishertestR.Pval)
    for ur = 1 : size(FishertestR.Pval{p},2)
        joint_idx = strcmp(FishertestR.rowdes.group,'Joint');
        pvals_corrected_joint = mafdr(FishertestR.Pval{p}(joint_idx,ur), 'BHFDR',true);
        count_joint(p,ur) = sum(pvals_corrected_joint<0.05);
        %datasets_joint(p,ur,:) = pvals_corrected_joint<0.05;
        
        muscle_idx = strcmp(FishertestR.rowdes.group,'Muscle');
        pvals_corrected_muscle = mafdr(FishertestR.Pval{p}(muscle_idx,ur), 'BHFDR',true);
        count_muscle(p,ur) = sum(pvals_corrected_muscle<0.05);
        %datasets_muscle(p,ur,:) = pvals_corrected_muscle<0.05;

    end
end

CombinedOverJointandMuscle.count_joint = count_joint;
CombinedOverJointandMuscle.count_muscle = count_muscle;

CombinedOverJointandMuscle
```


```text:Output
CombinedOverJointandMuscle = 
        cols_URs: {31x1 cell}
        rows_SPs: {12x1 cell}
     CombP_joint: [12x31 double]
    CombP_muscle: [12x31 double]
       CombP_all: [12x31 double]
       FDR_Joint: [12x31 double]
      FDR_Muscle: [12x31 double]
         FDR_All: [12x31 double]
     count_joint: [12x31 double]
    count_muscle: [12x31 double]

```

  
# Final output

```matlab:Code

CombinedPval  = table(repmat(CombinedOverJointandMuscle.rows_SPs,size(CombinedOverJointandMuscle.CombP_joint,2),1),...
                    repelem(CombinedOverJointandMuscle.cols_URs,size(CombinedOverJointandMuscle.CombP_joint,1)),...
                    reshape(CombinedOverJointandMuscle.CombP_joint,size(CombinedOverJointandMuscle.CombP_joint,1)*size(CombinedOverJointandMuscle.CombP_joint,2),1),...
                    reshape(CombinedOverJointandMuscle.CombP_muscle,size(CombinedOverJointandMuscle.CombP_muscle,1)*size(CombinedOverJointandMuscle.CombP_muscle,2),1),...
                    reshape(CombinedOverJointandMuscle.CombP_all,size(CombinedOverJointandMuscle.CombP_all,1)*size(CombinedOverJointandMuscle.CombP_all,2),1),...
                    reshape(CombinedOverJointandMuscle.FDR_Joint,size(CombinedOverJointandMuscle.FDR_Joint,1)*size(CombinedOverJointandMuscle.FDR_Joint,2),1),...
                    reshape(CombinedOverJointandMuscle.FDR_Muscle,size(CombinedOverJointandMuscle.FDR_Muscle,1)*size(CombinedOverJointandMuscle.FDR_Muscle,2),1),...
                    reshape(CombinedOverJointandMuscle.FDR_All,size(CombinedOverJointandMuscle.FDR_All,1)*size(CombinedOverJointandMuscle.FDR_All,2),1),...
                    reshape(CombinedOverJointandMuscle.count_joint,size(CombinedOverJointandMuscle.count_joint,1)*size(CombinedOverJointandMuscle.count_joint,2),1),...
                    reshape(CombinedOverJointandMuscle.count_muscle,size(CombinedOverJointandMuscle.count_muscle,1)*size(CombinedOverJointandMuscle.count_muscle,2),1),...
                    'variablenames',{'SP','UR','CombinedP_Joint','CombinedP_Muscle', 'CombinedP_All', 'FDR_Joint', 'FDR_Muscle', 'FDR_All', 'count_Joint', 'count_Muscle'});
head(CombinedPval, 5)
```

| |SP|UR|CombinedP_Joint|CombinedP_Muscle|CombinedP_All|FDR_Joint|FDR_Muscle|FDR_All|count_Joint|count_Muscle|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1|'1.1'|'APOE'|1.0000e+00|4.5845e-20|9.1843e-09|1.0000e+00|1.4212e-18|2.8471e-07|0|7|
|2|'1.10'|'APOE'|1.0000e+00|9.9784e-01|1.0000e+00|1.0000e+00|1.0000e+00|1.0000e+00|0|0|
|3|'1.2'|'APOE'|1.0000e+00|1.9419e-22|1.9182e-10|1.0000e+00|1.5050e-21|1.4866e-09|0|7|
|4|'1.3'|'APOE'|1.0000e+00|1.0711e-26|1.3463e-13|1.0000e+00|1.1068e-25|1.0434e-12|0|7|
|5|'1.4'|'APOE'|1.0000e+00|2.6564e-24|8.3587e-12|1.0000e+00|2.7449e-23|6.4780e-11|0|7|

  
### Save the data

\hfill \break


```matlab:Code
%writetable(CombinedPval,path_output)

```

