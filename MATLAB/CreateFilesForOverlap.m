function udis_condition_sumary = CreateFilesForOverlap(FN)


FN.GWAS_key = cellfun(@(x,y,z) sprintf('%s_%s_%s',x,y,z),FN.GSE,FN.disease_short,FN.tissue_detail,'UniformOutput',0);

FN.disease_short(strcmp(FN.GSE,'GSE81071') & strcmp(FN.disease_short,'DLE'))={'lupus'};
FN.disease_short(strcmp(FN.GSE,'GSE81071') & strcmp(FN.disease_short,'SCLE'))={'lupus'};
FN.disease_short(strcmp(FN.GSE,'GSE112943') & strcmp(FN.disease_short,'subacute cutaneous'))={'lupus'};
FN.disease_short(strcmp(FN.GSE,'GSE32591') & strcmp(FN.disease_short,'LN'))={'lupus'};
FN.disease_short(strcmp(FN.GSE,'GSE148810') & strcmp(FN.disease_short,'cSLE'))={'lupus'};
FN.disease_short(strcmp(FN.GSE,'GSE112943') & strcmp(FN.disease_short,'kidney'))={'lupus'};


FN.disease_key = FN.disease_short;
FN.disease_key(strcmp(FN.inflammed,'yes')) = cellfun(@(x) sprintf('%s_inflamed',x),FN.disease_key(strcmp(FN.inflammed,'yes')),'UniformOutput',0);
%FN.disease_key(~strcmp(FN.inflammed,'yes'));
FN.disease_key(strcmp(FN.inflammed,'no')) = cellfun(@(x) sprintf('%s_non-inflamed',x),FN.disease_key(strcmp(FN.inflammed,'no')),'UniformOutput',0);
%FN.disease_key(~strcmp(FN.inflammed,'no'));


udis_condition = table(unique(FN.disease_key),'variablenames',{'dis'});

for zz = 1 : length(udis_condition.dis)
	hv = FN.Pathways(strcmp(FN.disease_key,udis_condition.dis{zz}));
    if length(hv)>1
        for g = 1 : length(hv)
            if g==1
                hv{g}.logpvalue = 10.^-hv{g}.logpvalue;
                phv = hv{g}(:,{'IngenuityCanonicalPathways','logpvalue'});
            else
                hv{g}.logpvalue = 10.^-hv{g}.logpvalue;
                phv = outerjoin(phv,hv{g}(:,{'IngenuityCanonicalPathways','logpvalue'}),'key','IngenuityCanonicalPathways','MergeKeys',1);
            end
        end
        hv = phv(:,'IngenuityCanonicalPathways');
        phv = table2array(phv(:,2:end));

        % combine P-values:
        for hj = 1 : size(phv,1)
            phvhv = phv(hj,:);
            phvhv(isnan(phvhv))=[];
            chi = -2*sum(log(phvhv),2);
            hv.logpvalue(hj) = 1-chi2cdf(chi,2*size(phvhv,2));
        end
    else
        hv = hv{1};
        hv.logpvalue = 10.^-hv.logpvalue;
    end
    udis_condition.Paths{zz} = hv;
    clear phv hv
end
% save ('EnrichedPathwaysPerDisease_FDR_values_not_corrected','udis_condition')
for z=1:length(udis_condition.dis)
    if z == 1
        udis_condition_sumary = [table(repmat(udis_condition.dis(z),size(udis_condition.Paths{z},1),1),'VariableNames',{'disease'}),udis_condition.Paths{z}(:,{'IngenuityCanonicalPathways','logpvalue'})];
    else
        udis_condition_sumary = [udis_condition_sumary;...
            [table(repmat(udis_condition.dis(z),size(udis_condition.Paths{z},1),1),'VariableNames',{'disease'}),udis_condition.Paths{z}(:,{'IngenuityCanonicalPathways','logpvalue'})]];
    end
end
udis_condition_sumary.Properties.VariableNames(strcmp(udis_condition_sumary.Properties.VariableNames,'logpvalue')) = {'pvalue'};

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% uprograms = [{'P1'};{'P2'};unique(P.subprograms)];
% for z = 1 : length(uprograms)
%    if  strncmp(uprograms{z},'P',1)
%        paths_in_program = P.IngenuityCanonicalPathways(P.clustersGlobal == str2double(uprograms{z}(2)));
%    else
%        paths_in_program = P.IngenuityCanonicalPathways(strcmp(P.subprograms, uprograms{z}));
%    end
%    for zz = 1 : length(udis_condition.dis)
%       paths_in_disease = udis_condition.Paths{zz};
%      % background = P.IngenuityCanonicalPathways(ismember(P.IngenuityCanonicalPathways,paths_in_disease.IngenuityCanonicalPathways))
%       paths_in_disease(paths_in_disease.bhfdr>=0.05,:)=[];
%       paths_in_disease(~ismember(paths_in_disease.IngenuityCanonicalPathways,P.IngenuityCanonicalPathways),:)=[];
%       paths_in_disease = paths_in_disease.IngenuityCanonicalPathways;
%       
%       a = sum(ismember(paths_in_program,paths_in_disease));
%       b = sum(~ismember(paths_in_program,paths_in_disease));
%       c = sum(~ismember(paths_in_disease,paths_in_program));
%       d = sum(~ismember(P.IngenuityCanonicalPathways,[paths_in_program;paths_in_disease]));
%       [~, p(z,zz)] = fishertest([a b; c d],'tail','right');
%       clear a b c d
%    end
% end
% 
% writetable(array2table(p,'variablenames',udis_condition.dis,'rownames',uprograms),...
%     sprintf('%s%s.txt','/Users/danga10/Documents/CIA/DS_UR_Program_test/tree_structures_code_v2/overlap_diseases_programs/IMIDs_pathway_overlap_with_SPs_',selected_disease))
% 
% 
% p2 = reshape(p,1, size(p,1)*size(p,2));
% % p2(1:length(p(:,1)))'==p(:,1)
% % p2(1:25)'==p(:,1)
% % p2(26:50)'==p(:,2)
% % p2(26:50)'==p(:,3)
% if ~strcmp(selected_disease,'ALL') & ~strcmp(selected_disease,'tissue')
%     p2 = table(repmat(uprograms,size(p,2),1), reshape(repmat(udis_condition.dis',size(p,1),1),1,size(p,1)*size(p,2))',...
%         reshape(repmat(udis_condition.GWAS_key',size(p,1),1),1,size(p,1)*size(p,2))', p2','variablenames',{'program','udis_condition','GWAS_key','p_val'});
% else
%     p2 = table(repmat(uprograms,size(p,2),1), reshape(repmat(udis_condition.dis',size(p,1),1),1,size(p,1)*size(p,2))',p2','variablenames',{'program','udis_condition','p_val'});
% end
% % test_SP = '2.1'
% % testdua = p2(strcmp(p2.program,test_SP),:);
% % sum(~strcmp(testdua.udis_condition,udis_condition.dis))
% % sum(testdua.p_val~=p(strcmp(uprograms,test_SP),:)')
% p2.bhfdr=mafdr(p2.p_val);
% writetable(p2,...
%     sprintf('%s%s.txt','/Users/danga10/Documents/CIA/DS_UR_Program_test/tree_structures_code_v2/overlap_diseases_programs/IMIDs_pathway_overlap_with_SPs_reshaped_for_dot_plot',selected_disease))
