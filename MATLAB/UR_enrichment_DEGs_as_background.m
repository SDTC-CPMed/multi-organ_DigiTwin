clear all
close all
clc

% scRNAseq:
savename = 'CIA';
% load DEGs
PATHH = '../data/CIA_mouse_DEGs/';
Fn = struct2cell(dir(PATHH));
Fn(2:end,:)=[];
Fn(~contains(Fn,'.csv'))=[];
Fn = table(Fn','variablenames',{'fn'});
Fn.Cell = cellfun(@(x) x{5}, regexp(Fn.fn,'\_','split'),'UniformOutput',0);
Fn.Organ = cellfun(@(x) x{6}(1:end-4), regexp(Fn.fn,'\_','split'),'UniformOutput',0);
Fn(~ismember(Fn.Organ,{'Joint','Muscle'}),:)=[];

translation = readtable('../data/Connective_pathway_analysis/orthologous_translation_NicheNet_analysis.txt');
for z = 1 : length(Fn.fn)
   hv =readtable(sprintf('%s%s',PATHH,Fn.fn{z}));
   hv(~strcmp(hv.is_de_fdr_0_05,'True'),:)=[];
   xx= unique(translation(ismember(translation.Mouse_gene_name,hv.Var1),:));
   %sprintf('/n no of DEGs with no orthologs %s; total number of DEGs: %s',num2str(sum(~ismember(hv.Var1,xx.Mouse_gene_name))),num2str(length(hv.Var1)));
   Fn.DEG{z} = xx.Gene_name;
   clear hv xx
end
Fn.cell_dataset = cellfun(@(x) x{5},regexp(Fn.fn,'\_','split'),'UniformOutput',0);
Fn.Properties.VariableNames(strcmp(Fn.Properties.VariableNames,'Organ'))={'Organ_condition'};
Fn.group = repmat({'Inflamed'},length(Fn.fn),1);
Fn.group(strcmp(Fn.Organ_condition,'Muscle'))={'Noninflamed'};
Fn.key = cellfun(@(x,y) sprintf('%s***%s',x,y),Fn.cell_dataset,Fn.group,'UniformOutput',0);

Fn = Fn(:,{'key','cell_dataset','group','DEG'});
%% load SPs %% add subprograms of program 1!!!!
SPhv = readtable('../data/Connective_pathway_analysis/TreeStructure_nodes2_scRNAseq.txt');
SPhv.subprograms = cellfun(@(x) sprintf('1.%s',x), strtrim(cellstr(num2str(SPhv.subclusters))),'UniformOutput',0);
SP = table(unique(SPhv.subprograms),'variablenames',{'SP'});
for z = 1 : length(SP.SP)
    SP.AllMolecules{z} = unique(regexp(strjoin(SPhv.AllMolecules(strcmp(SPhv.subprograms,SP.SP{z})),','),'\,','split'));
end

SP.SP(length(SP.SP)+1)={'P1'};
SP.AllMolecules{end} = unique(regexp(strjoin(SPhv.AllMolecules(SPhv.clustersGlobal==1),','),'\,','split'));
SP.SP(length(SP.SP)+1)={'P2'};
SP.AllMolecules{end} = unique(regexp(strjoin(SPhv.AllMolecules(SPhv.clustersGlobal==2),','),'\,','split'));

%% load URs: "Just remember, the file needs to be sorted for: only positive PCCs and only intra-organ interactions"
UR = readtable('../data/Connective_pathway_analysis/all_outgoing_ligand_activity.txt');
UR(UR.pearson<0,:)=[];
UR.senderOrgan = cellfun(@(x) x{2}, regexp(UR.Sender,'\_','split'),'UniformOutput',0);
UR.targetOrgan = cellfun(@(x) x{2}, regexp(UR.Target,'\_','split'),'UniformOutput',0);
UR(~strcmp(UR.senderOrgan,UR.targetOrgan),:)=[];
UR(~ismember(UR.senderOrgan,{'Muscle','Joint'}),:)=[];
UR.targetcell = cellfun(@(x) x{1}, regexp(UR.Target,'\_','split'),'UniformOutput',0);
UR = unique(UR(:,{'test_ligand','targetcell','targetOrgan','target'}));
UR.Properties.VariableNames = {'UR','cell_dataset','condition','DS'};
UR.DS = regexp(UR.DS,'\/','split');
UR.group(strcmp(UR.condition,'Muscle')) = {'Noninflamed'};
UR.group(strcmp(UR.condition,'Joint')) = {'Inflamed'};
UR.key = cellfun(@(x,y) sprintf('%s***%s',x,y),UR.cell_dataset,UR.group,'UniformOutput',0);
UR = UR(:,{'UR','key','cell_dataset','group','DS'});


% UR_new = table( unique(UR.key),'variablenames',{'key'});
% for h = 1 : length(UR_new.key)
%     UR_new.DS{h} = unique(regexp(strjoin(UR.DS(strcmp(UR.key,UR_new.key{h})),'/'),'\/','split'));
% end


% %% IMIDs:
% savename = 'IMIDs';
% % load DEGs
% load('../data/Connective_pathway_analysis/AllDatasets_DEGs_martin.mat','AllDatasets')
% AllDatasets.key2 = cellfun(@(x,y,z) sprintf('%s_%s_%s',x,y,z),AllDatasets.gse,AllDatasets.disease_short,AllDatasets.tissue,'UniformOutput',0);
% AllDatasets(strcmp(AllDatasets.active_USE,'unknown'),:)=[];
% AllDatasets(strcmp(AllDatasets.key,'GSE40568_SS_labial_salivary_gland'),:)=[];
% 
% 
% load ../data/Connective_pathway_analysis/gene_info_type_of_gene_2020_03_03.mat
% gene_info = unique(gene_info(:,{'GeneID','Symbol','Synonyms','type_of_gene'}));
% datasets_group = {'GSE32924_uninflamed 12 vs HC 8_AD.csv','GSE16161_Skin AD 9 vs contol 9.csv','GSE32924_Skin_AD 13 VS_control8.csv','GSE16879_colon_CD.csv',...
%     'GSE179285_ascending descending colon_CD.csv','GSE179285_unflamed ascending descending colon 72 vs control 12_CD.csv','GSE16879_ileum_cd.csv',...
%     'GSE179285_Terminal ileum_CD.csv','GSE81071_DLE_vs_control.csv','GSE14905_psoriasis_non lesion skin 28 vs control 21.csv',...
%     'GSE1919_rheumatoid arthritis _vs healthy (5 vs 5).csv','GSE55235_rheumatoid arthritis _vs healthy (10 vs 10).csv','GSE81071_SCLE_vs_control.csv',...
%     'GSE40568_sjogren_vs_controllabial salivary glands.csv','GSE81292_SSc_ILD_lung.csv','GSE95065_SSC_skin.csv','GSE11223_descending_colon.csv',...
%     'GSE11223_sigmoid colon_UC.csv','GSE179285_sigmoid colon_UC.csv','GSE11223_Uninflamed 66 vs HC 69_UC.csv','GSE179285_inactive 32 vs contol 31_UC.csv',...
%     'GSE148810_childhood_onset_lupus_cSLE_skin.csv','GSE112943_kidney_lupus.csv','GSE112943_subacute cutaneous lupus.csv'};
% for z = 1 : length(AllDatasets.deg_file_name)
%     hv = readtable(join(['../data/AllDEGfilesMovedToOneFolder/', AllDatasets.deg_file_name{z}]));
%     if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE81071_DLE_vs_control.csv','GSE81071_SCLE_vs_control.csv','GSE95065_SSC_skin.csv'}))==1
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ENTREZ_GENE_ID'))={'Gene_ID'};
%     end
%     if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE148810_juvenile myositis_skin_1.csv','GSE148810_Nonlesional skin 6 vs HC 8_JM.csv'}))==1
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ORF'))={'Gene_symbol'};
%     end
%     if strcmp(AllDatasets.deg_file_name{z},'GSE32591_glomer_vs_contol_LN.csv')
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'Gene_Symbol'))={'Gene_symbol'};
%     end
%     if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE181318_skin_psoriatic 3 vs control3.csv','GSE66413_Pancreatic lymph nodes 13_ T1D vs healthy 3.csv'}))
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'GENE_SYMBOL'))={'Gene_symbol'};
%     end
%     if strcmp(AllDatasets.deg_file_name{z},'GSE176510_Sj”gren syndrome_keratocojunctivitis scicca.csv')
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ID'))={'Gene_symbol'};
%     end
%     if strcmp(AllDatasets.deg_file_name{z},'GSE81292_SSc_ILD_lung.csv')
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ORF'))={'Gene_ID'};
%     end
%     if strcmp(AllDatasets.deg_file_name{z},'GSE148810_childhood_onset_lupus_cSLE_skin.csv')
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ORF'))={'Gene_symbol'};
%         hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'SPOT_ID'))={'Gene_ID'};
%     end 
%     if sum(strcmp(datasets_group,AllDatasets.deg_file_name{z}))==1
%         genesfound = unique(hv.Gene_ID(~isnan(hv.Gene_ID)));
%         genesfound(~ismember(genesfound,gene_info.GeneID),:)=[];
%         AllDatasets.measured_genes{z}  = unique(gene_info.Symbol(ismember(gene_info.GeneID,genesfound)));
%         hv(hv.adj_P_Val>=0.05,:)=[];
%         genesfound(~ismember(genesfound,hv.Gene_ID),:)=[];
%         AllDatasets.DEG_genes{z}  = unique(gene_info.Symbol(ismember(gene_info.GeneID,genesfound)));
%         clear genesfound
%     elseif sum(strcmp(AllDatasets.deg_file_name{z},{'GSE75214_inactive_vs_normal_CD_16_11_ileum.csv','GSE148810_juvenile myositis_skin_1.csv',...
%             'GSE148810_Nonlesional skin 6 vs HC 8_JM.csv','GSE32591_glomer_vs_contol_LN.csv','GSE181318_skin_psoriatic 3 vs control3.csv',...
%             'GSE176510_Sj”gren syndrome_keratocojunctivitis scicca.csv','GSE75214_inactive_vs_normal_23_11_UC_colon.csv','GSE66413_Pancreatic lymph nodes 13_ T1D vs healthy 3.csv'}))==1
%         genesfound = unique(hv.Gene_symbol(~strcmp(hv.Gene_symbol,'')));
%         %%% For cases of genes with alternative names 
%         for i = 1:length(genesfound)
%             genes = regexp(genesfound(i),'///','split');
%             for j = 1:length(genes{1,1})
%                 if ismember(genes{1,1}(j), gene_info.Symbol)
%                     genesfound(i) = genes{1,1}(j);
%                 end
%             end
%         end       
%         genesfound(~ismember(genesfound,gene_info.Symbol),:)=[]; 
%         AllDatasets.measured_genes{z} = genesfound;
%         hv(hv.adj_P_Val>=0.05,:)=[];
%         genesfound = unique(hv.Gene_symbol(~strcmp(hv.Gene_symbol,'')));
%         %%% For cases of genes with alternative names
%         for i = 1:length(genesfound)
%             genes = regexp(genesfound(i),'///','split');
%             for j = 1:length(genes{1,1})
%                 if ismember(genes{1,1}(j), gene_info.Symbol)
%                     genesfound(i) = genes{1,1}(j);
%                 end
%             end
%         end
%         genesfound(~ismember(genesfound,gene_info.Symbol),:)=[];
%         AllDatasets.DEG_genes{z}  = unique(genesfound);
%         clear genesfound
%     end
% end
% 
% AllDatasets.broad_label = repmat({''},length(AllDatasets.key),1);
% AllDatasets.broad_label(strcmp(AllDatasets.active_USE,'yes')) = cellfun(@(x) sprintf('%s_active',x),AllDatasets.disease_short(strcmp(AllDatasets.active_USE,'yes')),'UniformOutput',0);
% AllDatasets(strcmp(AllDatasets.active_USE,'yes'),{'broad_label','disease_short'})
% AllDatasets.broad_label(strcmp(AllDatasets.active_USE,'no')) = cellfun(@(x) sprintf('%s_inactive',x),AllDatasets.disease_short(strcmp(AllDatasets.active_USE,'no')),'UniformOutput',0);
% AllDatasets(strcmp(AllDatasets.active_USE,'no'),{'broad_label','disease_short'})
% 
% AllDatasets.DEG=[];
% 
% AllDatasets.group = repmat({'Inflamed'},length(AllDatasets.DEG_genes),1);
% AllDatasets.group(strcmp(AllDatasets.active_USE,'no'))={'Noninflamed'};
% AllDatasets.Properties.VariableNames(strcmp(AllDatasets.Properties.VariableNames,'DEG_genes'))={'DEG'};
% AllDatasets.Properties.VariableNames(strcmp(AllDatasets.Properties.VariableNames,'disease_short'))={'cell_dataset'};
% AllDatasets.key = cellfun(@(x,y) sprintf('%s***%s',x,y),AllDatasets.file_name,AllDatasets.group,'UniformOutput',0);
% 
% broad_label = AllDatasets.broad_label;
% broad_label(10) = {'lupus_active'};
% broad_label(13) = {'lupus_active'};
% broad_label(18) = {'lupus_active'};
% broad_label(30) = {'lupus_active'};
% AllDatasets.broad_label = broad_label;
% 
% Fn = AllDatasets(:,{'key','file_name','group','DEG', 'broad_label'});
% clearvars -except Fn savename AllDatasets gene_info
% 
% 
% %% load SPs %% add subprograms of program 1
% P2 = readtable('../data/Connective_pathway_analysis/TreeStructure_nodes2_CLUSTER2_AID_noblood.txt');
% P2.subprograms = cellfun(@(x) sprintf('2.%s',x), strtrim(cellstr(num2str(P2.subclusters))),'UniformOutput',0);
% P1 = readtable('../data/Connective_pathway_analysis/TreeStructure_nodes2_AID_noblood.txt');
% P1.subprograms = cellfun(@(x) sprintf('1.%s',x), strtrim(cellstr(num2str(P1.subclusters))),'UniformOutput',0);
% P = [P1(:,{'clustersGlobal','subprograms','AllMolecules'});P2(:,{'clustersGlobal','subprograms','AllMolecules'})];
% Programs = table(unique(P.subprograms),'variablenames',{'Program'});
% for z = 1 : length(Programs.Program)
%     Programs.AllMolecules{z} = unique(regexp(strjoin(P.AllMolecules(strcmp(P.subprograms,Programs.Program{z})),','),'\,','split'));
% end
% 
% Programs.Program(length(Programs.Program)+1)={'P1'};
% Programs.AllMolecules{end} = unique(regexp(strjoin(P.AllMolecules(P.clustersGlobal==1),','),'\,','split'));
% 
% Programs.Program(length(Programs.Program)+1)={'P2'};
% Programs.AllMolecules{end} = unique(regexp(strjoin(P.AllMolecules(P.clustersGlobal==2),','),'\,','split'));
% 
% SP=Programs;
% SP.Properties.VariableNames(strcmp(SP.Properties.VariableNames,'Program'))={'SP'};
% clearvars -except SP Fn savename AllDatasets gene_info
% 
% %% load URs: 
% AllDatasets = AllDatasets(:,{'file_name','key','cell_dataset','tissue','gse','active_USE','UR', 'measured_genes','DEG','group'});
% 
% for z = 1 : length(AllDatasets.key)
%     AllDatasets.UR{z}(AllDatasets.UR{z}.p_valueofoverlap>0.05,:)=[]; % I am using all to get an overview on ALL possible target genes of an UR
%     %AllDatasets.UR{z}(~ismember(AllDatasets.UR{z}.upstreamregulator,AllDatasets.DEG{z}),:) =[]; %Check if it is also DEG
%     hv = [table(repmat(AllDatasets.file_name(z),length(AllDatasets.UR{z}.upstreamregulator),1),...
%         repmat(AllDatasets.group(z),length(AllDatasets.UR{z}.upstreamregulator),1),'variablenames',{'file_name','group'}),...
%         AllDatasets.UR{z}(:,{'upstreamregulator','targetmoleculesindataset'})];
%     if z == 1
%         allURhv = hv;
%     else
%         
%         allURhv = [allURhv;hv];
%     end
% end
% 
% allURhv.targetmoleculesindataset =  regexp(allURhv.targetmoleculesindataset,'\,','split');
% 
% %%% Check if DSs are among DEGs
% for i = 1:length(allURhv.targetmoleculesindataset)
%     DEGs = AllDatasets.DEG{strcmp(AllDatasets.file_name, allURhv.file_name(i))};
%     allURhv.targetmoleculesindataset{i,1}(~ismember(allURhv.targetmoleculesindataset{i,1},DEGs)) = [];
% end
%     
% allURhv.Properties.VariableNames(strcmp(allURhv.Properties.VariableNames,'targetmoleculesindataset'))={'DS'};
% allURhv.Properties.VariableNames(strcmp(allURhv.Properties.VariableNames,'upstreamregulator'))={'UR'};
% allURhv.key  = cellfun(@(x,y) sprintf('%s***%s',x,y),allURhv.file_name,allURhv.group,'UniformOutput',0);
% UR = allURhv(:,{'UR','key','file_name','group','DS'});
% clearvars -except SP Fn savename UR
% 


%% %%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%% %% %%%%%%%%%%%%%%%%%%%%% %%
%% test Fisher test enrichment of UR DS in SPs:
digitsOld = digits(100);
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

FishertestR.cellsdes = SP;
FishertestR.rowdes = Fn;
FishertestR.coldes_URs = uUR;
FishertestR.moreURinfo = UR;
FishertestR.Pval = FPval;
FishertestR.OR = FOR;


%% combine pvalues over inflamed and noninflamed:

for p = 1 : length(FishertestR.Pval)
    for ur = 1 : size(FishertestR.Pval{p},2)
        inflamed_idx = strcmp(FishertestR.rowdes.group,'Inflamed');
        chi = -2*sum(log(FishertestR.Pval{p}(inflamed_idx,ur)));
        CombinedPval_Inflamed(p,ur) = my_chi2cdf(chi,2*sum(inflamed_idx));
        
        noninflamed_idx = strcmp(FishertestR.rowdes.group,'Noninflamed');
        chi = -2*sum(log(FishertestR.Pval{p}(noninflamed_idx,ur)));
        CombinedPval_Noninflamed(p,ur) = my_chi2cdf(chi,2*sum(noninflamed_idx));
        
        chi = -2*sum(log(FishertestR.Pval{p}(:,ur)));
        CombinedPval_all(p,ur) = my_chi2cdf(chi,2*32);
    end
end


CombinedOverInfandNoninf.cols_URs = FishertestR.coldes_URs;
CombinedOverInfandNoninf.rows_SPs = FishertestR.cellsdes.SP;
CombinedOverInfandNoninf.CombP_Inf = CombinedPval_Inflamed;
CombinedOverInfandNoninf.CombP_Noninf = CombinedPval_Noninflamed;
CombinedOverInfandNoninf.CombP_all = CombinedPval_all;


%% FDR correction over Programs 
for i = 1:length(CombinedOverInfandNoninf.CombP_Inf(:,1))
    qval_Inf(i,:) = mafdr(CombinedOverInfandNoninf.CombP_Inf(i,:), 'BHFDR',true);
    qval_Noninf(i,:) = mafdr(CombinedOverInfandNoninf.CombP_Noninf(i,:), 'BHFDR',true);
    qval_All(i,:) = mafdr(CombinedOverInfandNoninf.CombP_all(i,:), 'BHFDR',true);

end

CombinedOverInfandNoninf.qval_Inf = qval_Inf;
CombinedOverInfandNoninf.qval_Noninf = qval_Noninf;
CombinedOverInfandNoninf.qval_All = qval_All;


CombinedOverInfandNoninf.cols_URs(CombinedOverInfandNoninf.CombP_Inf(strcmp(CombinedOverInfandNoninf.rows_SPs,'1.6'),:)<0.05)
CombinedOverInfandNoninf.cols_URs(CombinedOverInfandNoninf.CombP_Noninf(strcmp(CombinedOverInfandNoninf.rows_SPs,'1.6'),:)<0.05)

%% Restructure the data for each dataset
for p = 1 : length(FishertestR.Pval)
    for ur = 1 : size(FishertestR.Pval{p},2)
        unique_diseases = unique(FishertestR.rowdes.key);
        for disease = 1 : length(unique_diseases)
            disease_idx = strcmp(FishertestR.rowdes.key,unique_diseases(disease));
            CombinedPval_Dataset(p,ur,disease) = FishertestR.Pval{p}(disease_idx,ur);
            CombinedOR_dataset(p,ur,disease) = FishertestR.OR{p}(disease_idx,ur);
        end
    end
end

%% Combine pvals for each disease
if savename == "IMIDs"
    for p = 1 : length(FishertestR.Pval)
        for ur = 1 : size(FishertestR.Pval{p},2)
            unique_diseases = unique(FishertestR.rowdes.broad_label);
            for disease = 1 : length(unique_diseases)
                disease_idx = strcmp(FishertestR.rowdes.broad_label,unique_diseases(disease));
                chi = -2*sum(log(FishertestR.Pval{p}(disease_idx,ur)));
                CombinedPval_Disease(p,ur,disease) = my_chi2cdf(chi,2*sum(disease_idx));           
            end
        end
    end
end




%% count pvalues over inflamed and noninflamed:

% For each program and UR, do FDR correction over all datasets 
for p = 1 : length(FishertestR.Pval)
    for ur = 1 : size(FishertestR.Pval{p},2)
        inflamed_idx = strcmp(FishertestR.rowdes.group,'Inflamed');
        pvals_corrected_Inf = mafdr(FishertestR.Pval{p}(inflamed_idx,ur), 'BHFDR',true);
        count_Inf(p,ur) = sum(pvals_corrected_Inf<0.05);
        datasets_Inf(p,ur,:) = pvals_corrected_Inf<0.05;
        
        noninflamed_idx = strcmp(FishertestR.rowdes.group,'Noninflamed');
        pvals_corrected_Noninf = mafdr(FishertestR.Pval{p}(noninflamed_idx,ur), 'BHFDR',true);
        count_Noninf(p,ur) = sum(pvals_corrected_Noninf<0.05);
        datasets_Noninf(p,ur,:) = pvals_corrected_Noninf<0.05;

    end
end

CombinedOverInfandNoninf.count_Inf = count_Inf;
CombinedOverInfandNoninf.count_Noninf = count_Noninf;

if savename == "IMIDs"
    CombinedPval  = table(repmat(CombinedOverInfandNoninf.rows_SPs,size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    repelem(CombinedOverInfandNoninf.cols_URs,size(CombinedOverInfandNoninf.CombP_Inf,1)),...
                    reshape(CombinedOverInfandNoninf.CombP_Inf,size(CombinedOverInfandNoninf.CombP_Inf,1)*size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_Noninf,size(CombinedOverInfandNoninf.CombP_Noninf,1)*size(CombinedOverInfandNoninf.CombP_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_all,size(CombinedOverInfandNoninf.CombP_all,1)*size(CombinedOverInfandNoninf.CombP_all,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_Inf,size(CombinedOverInfandNoninf.qval_Inf,1)*size(CombinedOverInfandNoninf.qval_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_Noninf,size(CombinedOverInfandNoninf.qval_Noninf,1)*size(CombinedOverInfandNoninf.qval_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_All,size(CombinedOverInfandNoninf.qval_All,1)*size(CombinedOverInfandNoninf.qval_All,2),1),...
                    reshape(CombinedOverInfandNoninf.count_Inf,size(CombinedOverInfandNoninf.count_Inf,1)*size(CombinedOverInfandNoninf.count_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.count_Noninf,size(CombinedOverInfandNoninf.count_Noninf,1)*size(CombinedOverInfandNoninf.count_Noninf,2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,1)),size(squeeze(CombinedPval_Disease(:,:,1)),1)*size(squeeze(CombinedPval_Disease(:,:,1)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,2)),size(squeeze(CombinedPval_Disease(:,:,2)),1)*size(squeeze(CombinedPval_Disease(:,:,2)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,3)),size(squeeze(CombinedPval_Disease(:,:,3)),1)*size(squeeze(CombinedPval_Disease(:,:,3)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,4)),size(squeeze(CombinedPval_Disease(:,:,4)),1)*size(squeeze(CombinedPval_Disease(:,:,4)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,5)),size(squeeze(CombinedPval_Disease(:,:,5)),1)*size(squeeze(CombinedPval_Disease(:,:,5)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,6)),size(squeeze(CombinedPval_Disease(:,:,6)),1)*size(squeeze(CombinedPval_Disease(:,:,6)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,7)),size(squeeze(CombinedPval_Disease(:,:,7)),1)*size(squeeze(CombinedPval_Disease(:,:,7)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,8)),size(squeeze(CombinedPval_Disease(:,:,8)),1)*size(squeeze(CombinedPval_Disease(:,:,8)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,9)),size(squeeze(CombinedPval_Disease(:,:,9)),1)*size(squeeze(CombinedPval_Disease(:,:,9)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,10)),size(squeeze(CombinedPval_Disease(:,:,10)),1)*size(squeeze(CombinedPval_Disease(:,:,10)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,11)),size(squeeze(CombinedPval_Disease(:,:,11)),1)*size(squeeze(CombinedPval_Disease(:,:,11)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,12)),size(squeeze(CombinedPval_Disease(:,:,12)),1)*size(squeeze(CombinedPval_Disease(:,:,12)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,13)),size(squeeze(CombinedPval_Disease(:,:,13)),1)*size(squeeze(CombinedPval_Disease(:,:,13)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,14)),size(squeeze(CombinedPval_Disease(:,:,14)),1)*size(squeeze(CombinedPval_Disease(:,:,14)),2),1),...
                    reshape(squeeze(CombinedPval_Disease(:,:,15)),size(squeeze(CombinedPval_Disease(:,:,15)),1)*size(squeeze(CombinedPval_Disease(:,:,15)),2),1),...
                    'variablenames',{'SP','UR','CombinedP_Inflamed','CombinedP_Noninflamed', 'CombinedP_All', 'qval_Inflamed', 'qval_Noninflamed', 'qval_All', 'count_Inflamed', 'count_Noninflamed','AD_active','AD_inactive','CD_active','CD_inactive','JM_active','JM_inactive','PSO_active','PSO_inactive','RA_active','SS_active','SSc_active','UC_active','UC_inactive','at_risk_T1D_inactive','lupus_active'});

    CombinedPval_lupus  = table(repmat(CombinedOverInfandNoninf.rows_SPs,size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    repelem(CombinedOverInfandNoninf.cols_URs,size(CombinedOverInfandNoninf.CombP_Inf,1)),...
                    reshape(CombinedOverInfandNoninf.CombP_Inf,size(CombinedOverInfandNoninf.CombP_Inf,1)*size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_Noninf,size(CombinedOverInfandNoninf.CombP_Noninf,1)*size(CombinedOverInfandNoninf.CombP_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_Inf,size(CombinedOverInfandNoninf.qval_Inf,1)*size(CombinedOverInfandNoninf.qval_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_Noninf,size(CombinedOverInfandNoninf.qval_Noninf,1)*size(CombinedOverInfandNoninf.qval_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.count_Inf,size(CombinedOverInfandNoninf.count_Inf,1)*size(CombinedOverInfandNoninf.count_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.count_Noninf,size(CombinedOverInfandNoninf.count_Noninf,1)*size(CombinedOverInfandNoninf.count_Noninf,2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,10)),size(squeeze(CombinedPval_Dataset(:,:,10)),1)*size(squeeze(CombinedPval_Dataset(:,:,10)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,13)),size(squeeze(CombinedPval_Dataset(:,:,13)),1)*size(squeeze(CombinedPval_Dataset(:,:,13)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,18)),size(squeeze(CombinedPval_Dataset(:,:,18)),1)*size(squeeze(CombinedPval_Dataset(:,:,18)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,19)),size(squeeze(CombinedPval_Dataset(:,:,19)),1)*size(squeeze(CombinedPval_Dataset(:,:,19)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,20)),size(squeeze(CombinedPval_Dataset(:,:,20)),1)*size(squeeze(CombinedPval_Dataset(:,:,20)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,30)),size(squeeze(CombinedPval_Dataset(:,:,30)),1)*size(squeeze(CombinedPval_Dataset(:,:,30)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,31)),size(squeeze(CombinedPval_Dataset(:,:,31)),1)*size(squeeze(CombinedPval_Dataset(:,:,31)),2),1),...
                    reshape(squeeze(CombinedPval_Dataset(:,:,32)),size(squeeze(CombinedPval_Dataset(:,:,32)),1)*size(squeeze(CombinedPval_Dataset(:,:,32)),2),1),...
                   'variablenames',{'SP','UR','CombinedP_Inflamed','CombinedP_Noninflamed', 'qval_Inflamed', 'qval_Noninflamed', 'count_Inflamed', 'count_Noninflamed','URs_GSE81071_DLE_skin','URs_GSE32591_glomer_LN','URs_GSE81071_SCLE_skin','UR_GSE176510_Sjogren syndrome_keratocojunctivitis','URs_GSE40568_SJogren_labial salivary gland','UR_GSE148810_childhood-onset lupus(cSLE)_skin','UR_GSE112943_lupus_kidney','UR_GSE112943_lupus_subacute cutaneous'});

    %These are used in the Python code to construct heatmaps
    datasets_Inf = squeeze(datasets_Inf(7,:,:));
    datasets_Noninf = squeeze(datasets_Noninf(7,:,:));
    
    %writetable(CombinedPval_lupus,sprintf('../data/UR_analysis/UR_predictions_%s_lupus_SS',savename))
    %writematrix(datasets_Inf,sprintf('../data/UR_analysis/Datasets_Inf',savename))
    %writematrix(datasets_Noninf,sprintf('../data/UR_analysis/Datasets_Noninf',savename))
    %writecell(CombinedOverInfandNoninf.cols_URs,sprintf('../data/UR_analysis/URs',savename))
               
else 
    CombinedPval  = table(repmat(CombinedOverInfandNoninf.rows_SPs,size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    repelem(CombinedOverInfandNoninf.cols_URs,size(CombinedOverInfandNoninf.CombP_Inf,1)),...
                    reshape(CombinedOverInfandNoninf.CombP_Inf,size(CombinedOverInfandNoninf.CombP_Inf,1)*size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_Noninf,size(CombinedOverInfandNoninf.CombP_Noninf,1)*size(CombinedOverInfandNoninf.CombP_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_all,size(CombinedOverInfandNoninf.CombP_all,1)*size(CombinedOverInfandNoninf.CombP_all,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_Inf,size(CombinedOverInfandNoninf.qval_Inf,1)*size(CombinedOverInfandNoninf.qval_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_Noninf,size(CombinedOverInfandNoninf.qval_Noninf,1)*size(CombinedOverInfandNoninf.qval_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.qval_All,size(CombinedOverInfandNoninf.qval_All,1)*size(CombinedOverInfandNoninf.qval_All,2),1),...
                    reshape(CombinedOverInfandNoninf.count_Inf,size(CombinedOverInfandNoninf.count_Inf,1)*size(CombinedOverInfandNoninf.count_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.count_Noninf,size(CombinedOverInfandNoninf.count_Noninf,1)*size(CombinedOverInfandNoninf.count_Noninf,2),1),...
                    'variablenames',{'SP','UR','CombinedP_Joint','CombinedP_Muscle', 'CombinedP_All', 'qval_Joint', 'qval_Muscle', 'qval_All', 'count_Joint', 'count_Muscle'});

end
        
% CombinedOR  = table(repmat(CombinedOverInfandNoninf.rows_SPs,size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
%                 repelem(CombinedOverInfandNoninf.cols_URs,size(CombinedOverInfandNoninf.CombP_Inf,1)),...
%                 reshape(CombinedOverInfandNoninf.CombP_Inf,size(CombinedOverInfandNoninf.CombP_Inf,1)*size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
%                 reshape(CombinedOverInfandNoninf.CombP_Noninf,size(CombinedOverInfandNoninf.CombP_Noninf,1)*size(CombinedOverInfandNoninf.CombP_Noninf,2),1),...
%                 reshape(CombinedOverInfandNoninf.qval_Inf,size(CombinedOverInfandNoninf.qval_Inf,1)*size(CombinedOverInfandNoninf.qval_Inf,2),1),...
%                 reshape(CombinedOverInfandNoninf.qval_Noninf,size(CombinedOverInfandNoninf.qval_Noninf,1)*size(CombinedOverInfandNoninf.qval_Noninf,2),1),...
%                 reshape(CombinedOverInfandNoninf.count_Inf,size(CombinedOverInfandNoninf.count_Inf,1)*size(CombinedOverInfandNoninf.count_Inf,2),1),...
%                 reshape(CombinedOverInfandNoninf.count_Noninf,size(CombinedOverInfandNoninf.count_Noninf,1)*size(CombinedOverInfandNoninf.count_Noninf,2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,1)),size(squeeze(CombinedOR_Dataset(:,:,1)),1)*size(squeeze(CombinedOR_Dataset(:,:,1)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,2)),size(squeeze(CombinedOR_Dataset(:,:,2)),1)*size(squeeze(CombinedOR_Dataset(:,:,2)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,3)),size(squeeze(CombinedOR_Dataset(:,:,3)),1)*size(squeeze(CombinedOR_Dataset(:,:,3)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,4)),size(squeeze(CombinedOR_Dataset(:,:,4)),1)*size(squeeze(CombinedOR_Dataset(:,:,4)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,5)),size(squeeze(CombinedOR_Dataset(:,:,5)),1)*size(squeeze(CombinedOR_Dataset(:,:,5)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,6)),size(squeeze(CombinedOR_Dataset(:,:,6)),1)*size(squeeze(CombinedOR_Dataset(:,:,6)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,7)),size(squeeze(CombinedOR_Dataset(:,:,7)),1)*size(squeeze(CombinedOR_Dataset(:,:,7)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,8)),size(squeeze(CombinedOR_Dataset(:,:,8)),1)*size(squeeze(CombinedOR_Dataset(:,:,8)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,9)),size(squeeze(CombinedOR_Dataset(:,:,9)),1)*size(squeeze(CombinedOR_Dataset(:,:,9)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,10)),size(squeeze(CombinedOR_Dataset(:,:,10)),1)*size(squeeze(CombinedOR_Dataset(:,:,10)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,11)),size(squeeze(CombinedOR_Dataset(:,:,11)),1)*size(squeeze(CombinedOR_Dataset(:,:,11)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,12)),size(squeeze(CombinedOR_Dataset(:,:,12)),1)*size(squeeze(CombinedOR_Dataset(:,:,12)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,13)),size(squeeze(CombinedOR_Dataset(:,:,13)),1)*size(squeeze(CombinedOR_Dataset(:,:,13)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,14)),size(squeeze(CombinedOR_Dataset(:,:,14)),1)*size(squeeze(CombinedOR_Dataset(:,:,14)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,15)),size(squeeze(CombinedOR_Dataset(:,:,15)),1)*size(squeeze(CombinedOR_Dataset(:,:,15)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,16)),size(squeeze(CombinedOR_Dataset(:,:,16)),1)*size(squeeze(CombinedOR_Dataset(:,:,16)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,17)),size(squeeze(CombinedOR_Dataset(:,:,17)),1)*size(squeeze(CombinedOR_Dataset(:,:,17)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,18)),size(squeeze(CombinedOR_Dataset(:,:,18)),1)*size(squeeze(CombinedOR_Dataset(:,:,18)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,19)),size(squeeze(CombinedOR_Dataset(:,:,19)),1)*size(squeeze(CombinedOR_Dataset(:,:,19)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,20)),size(squeeze(CombinedOR_Dataset(:,:,20)),1)*size(squeeze(CombinedOR_Dataset(:,:,20)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,21)),size(squeeze(CombinedOR_Dataset(:,:,21)),1)*size(squeeze(CombinedOR_Dataset(:,:,21)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,22)),size(squeeze(CombinedOR_Dataset(:,:,22)),1)*size(squeeze(CombinedOR_Dataset(:,:,22)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,23)),size(squeeze(CombinedOR_Dataset(:,:,23)),1)*size(squeeze(CombinedOR_Dataset(:,:,23)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,24)),size(squeeze(CombinedOR_Dataset(:,:,24)),1)*size(squeeze(CombinedOR_Dataset(:,:,24)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,25)),size(squeeze(CombinedOR_Dataset(:,:,25)),1)*size(squeeze(CombinedOR_Dataset(:,:,25)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,26)),size(squeeze(CombinedOR_Dataset(:,:,26)),1)*size(squeeze(CombinedOR_Dataset(:,:,26)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,27)),size(squeeze(CombinedOR_Dataset(:,:,27)),1)*size(squeeze(CombinedOR_Dataset(:,:,27)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,28)),size(squeeze(CombinedOR_Dataset(:,:,28)),1)*size(squeeze(CombinedOR_Dataset(:,:,28)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,29)),size(squeeze(CombinedOR_Dataset(:,:,29)),1)*size(squeeze(CombinedOR_Dataset(:,:,29)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,30)),size(squeeze(CombinedOR_Dataset(:,:,30)),1)*size(squeeze(CombinedOR_Dataset(:,:,30)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,31)),size(squeeze(CombinedOR_Dataset(:,:,31)),1)*size(squeeze(CombinedOR_Dataset(:,:,31)),2),1),...
%                 reshape(squeeze(CombinedOR_Dataset(:,:,32)),size(squeeze(CombinedOR_Dataset(:,:,32)),1)*size(squeeze(CombinedOR_Dataset(:,:,32)),2),1),...
%                'variablenames',{'SP','UR','CombinedP_Inflamed','CombinedP_Noninflamed', 'qval_Inflamed', 'qval_Noninflamed', 'count_Inflamed', 'count_Noninflamed','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32'});

%writetable(CombinedOR,sprintf('../data/UR_analysis/UR_OR_datasets',savename))


writetable(CombinedPval,sprintf('../data/UR_analysis/UR_predictions_%s_disease_Pvals_all.xlsx',savename))



%% Get z_scores after section 1 is run

%for z = 1 : length(AllDatasets.key)
%
%    AllDatasets.UR{z}(AllDatasets.UR{z}.p_valueofoverlap>0.05,:)=[]; % I am using all to get an overview on ALL possible target genes of an UR
%    %AllDatasets.UR{z}(~ismember(AllDatasets.UR{z}.upstreamregulator,AllDatasets.DEG{z}),:) =[];
%    z_scores = AllDatasets.UR{z}(:,1:2);
%    %writetable(z_scores,sprintf(join(['../data/UR_analysis/z_scores/',string(AllDatasets.deg_file_name{z})]),savename))
%        
%        
%end