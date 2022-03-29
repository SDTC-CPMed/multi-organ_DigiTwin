
%% read data GSE92415
%URs = readtable('GSE92415/nonresponders_URs');
%DEGs = readtable('GSE92415/GSE92415.untreatedNonResponder_vs_control_all_significant.csv');
URs = readtable('../data/GSE92415/responders_URs');
DEGs = readtable('../data/GSE92415/GSE92415.untreatedResponder_vs_control_all_significant.csv');

%% read data GSE73661
%URs = readtable('GSE73661/nonresponders_URs');
%DEGs = readtable('GSE73661/GSE73661.IFX_untreatedNonResponder_vs_control_all_significant.csv');
%URs = readtable('GSE73661/responders_URs');
%DEGs = readtable('GSE73661/GSE73661.IFX_untreatedResponder_vs_control_all_significant.csv');




%% preprocess
DEGs = transpose(DEGs.Gene_symbol);

URs = URs(strcmp(URs.Molecule_Type, "G-protein coupled receptor") |...
    strcmp(URs.Molecule_Type, 'cytokine') |...
    strcmp(URs.Molecule_Type, 'growth factor') |...
    strcmp(URs.Molecule_Type, 'ligand-dependent nuclear receptor') |...
    strcmp(URs.Molecule_Type, 'transmembrane receptor'),:);


%% load SPs %% add subprograms of program 1
P2 = readtable('TreeStructure_nodes2_CLUSTER2_AID_noblood.txt');
P2.subprograms = cellfun(@(x) sprintf('2.%s',x), strtrim(cellstr(num2str(P2.subclusters))),'UniformOutput',0);
P1 = readtable('TreeStructure_nodes2_AID_noblood.txt');
P1.subprograms = cellfun(@(x) sprintf('1.%s',x), strtrim(cellstr(num2str(P1.subclusters))),'UniformOutput',0);
P = [P1(:,{'clustersGlobal','subprograms','AllMolecules'});P2(:,{'clustersGlobal','subprograms','AllMolecules'})];
Programs = table(unique(P.subprograms),'variablenames',{'Program'});
for z = 1 : length(Programs.Program)
    Programs.AllMolecules{z} = unique(regexp(strjoin(P.AllMolecules(strcmp(P.subprograms,Programs.Program{z})),','),'\,','split'));
end

Programs.Program(length(Programs.Program)+1)={'P1'};
Programs.AllMolecules{end} = unique(regexp(strjoin(P.AllMolecules(P.clustersGlobal==1),','),'\,','split'));

Programs.Program(length(Programs.Program)+1)={'P2'};
Programs.AllMolecules{end} = unique(regexp(strjoin(P.AllMolecules(P.clustersGlobal==2),','),'\,','split'));

SP=Programs;
SP.Properties.VariableNames(strcmp(SP.Properties.VariableNames,'Program'))={'SP'};


%% test Fisher test enrichment of UR DS in SPs:
digitsOld = digits(100);
FPval = cell(length(SP.SP),1);
FOR = cell(length(SP.SP),1);
spgenes_list = strings(length(SP.SP),length(DEGs));

uUR = URs.Upstream_Regulator;
for p = 1 : length(SP.SP)
                 
       spgenes = SP.AllMolecules{p}(ismember(SP.AllMolecules{p},DEGs));
       for ur = 1 : length(uUR)          
                urds = URs.Target_Molecules_in_Dataset{ur}; 
                urds = regexp(urds,'\,','split');
                a = sum(ismember(spgenes,urds));
                b = sum(~ismember(spgenes,urds));
                c = sum(~ismember(urds,spgenes));
                d = sum(~ismember(DEGs,[spgenes,urds]));

                [~, FPval{p}(ur)] = fishertest([a b; c d],'tail','right');

                FOR{p}(ur) = (a*d)/(b*c);                
                clear urds a b c d
       end
       spgenes_list(p,1:length(spgenes)) = spgenes; 
end

FishertestR.cellsdes = SP;
FishertestR.coldes_URs = uUR;
FishertestR.moreURinfo = URs;
FishertestR.Pval = FPval;
FishertestR.OR = FOR;


%% FDR correction over Programs 
for i = 1:length(FishertestR.Pval(:,1))
    qval(i,:) = mafdr(FishertestR.Pval{i,:}, 'BHFDR',true);
end

sTable = array2table(qval,'RowNames',SP.SP,'VariableNames',uUR);

%writetable(sTable,sprintf('GSE73661/UR_predictions_nonresponders'), 'WriteRowNames',true)


%% DEGs for SPs
%this is only relevant for non-responders as we need SP specific DEGs only
%for non-respondes

sTable = array2table(spgenes_list,'RowNames',SP.SP);

%writetable(sTable,sprintf('GSE73661/DEGs_SP'), 'WriteRowNames',true)



