function [AllEnrichedPaths,filesForOverlap] = ReadInIPAPathwayEnrichments(FN,PATHH1)

cell_orgainzation = sort(FN.key(strcmp(FN.inflammed,'yes')));
cell_orgainzation_savename = cellfun(@(x) sprintf('%s_likeJoint',x),cell_orgainzation,'UniformOutput',0);
hv  = flip(sort(FN.key(strcmp(FN.inflammed,'no'))));
hv = cellfun(@(x) sprintf('%s_likeMuscle',x),hv,'UniformOutput',0);
cell_orgainzation =[cell_orgainzation;flip(sort(FN.key(strcmp(FN.inflammed,'no'))))];
cell_orgainzation_savename = [cell_orgainzation_savename;hv];

count=0;
for z = 1 : length(FN.file)
   hv = readtable(sprintf('%s/%s',PATHH1,FN.file{z}));
   if strcmp(hv.Properties.VariableNames{2},'x_log_p_value_')
      hv.Properties.VariableNames{2} ='logpvalue';
   end
   if strcmp(hv.Properties.VariableNames{4},'z_score')
      hv.Properties.VariableNames{4} ='zscore';
   end
   if strcmp(hv.Properties.VariableNames{1},'x_2000_2021QIAGEN_AllRightsReserved_')
       hv.Properties.VariableNames = strrep(strrep(strrep(strrep(table2array(hv(1,:)),'-',''),' ',''),'(',''),')','');
       hv(1,:)=[];
       if  contains(hv.logpvalue{3},'×')
           hvd = str2double(strrep(strrep(cellfun(@(x) x{1},regexp(hv.logpvalue,'\×','split'),'UniformOutput',0),',','.'),'−','-'));
           hvd2 = (cellfun(@(x) x{2},regexp(hv.logpvalue,'\^','split'),'UniformOutput',0));
           hvd2 = str2double(strrep(hvd2,'−','-'));
           hv.logpvalue = hvd.* 10.^hvd2;
           clear hvd hvd2
           
           hv.zscore = str2double(strrep(strrep(hv.zscore,',','.'),'−','-'));
           count=count+1;
       else           
           hv.logpvalue = str2double(strrep(hv.logpvalue,',','.'));
           hv.zscore = str2double(strrep(hv.zscore,',','.'));
       end
   elseif ~strcmp(hv.Properties.VariableNames{1},'IngenuityCanonicalPathways')
       error('Check headers!')
   end
   if iscell(hv.logpvalue)
       hv.logpvalue = str2double(strrep(hv.logpvalue,',','.'));
   end
   if iscell(hv.zscore)
       hv.zscore = str2double(strrep(hv.zscore,',','.'));
   end

%    hv(hv.logpvalue<-log10(0.05),:)=[]; % remove paths not enriched
   hv.IngenuityCanonicalPathways = strrep(hv.IngenuityCanonicalPathways,'Î±','α');
   hv.IngenuityCanonicalPathways = strrep(hv.IngenuityCanonicalPathways,'Î²','β');
   hv.IngenuityCanonicalPathways = strrep(hv.IngenuityCanonicalPathways,'Î³','γ');
   hv.IngenuityCanonicalPathways = strrep(hv.IngenuityCanonicalPathways,'Îº','κ');
   hv.IngenuityCanonicalPathways = strrep(hv.IngenuityCanonicalPathways,'Ïƒ','σ');

%    if z == 1 
%        AllEnrichedPaths = unique(hv.IngenuityCanonicalPathways);
%    else
%        AllEnrichedPaths = [AllEnrichedPaths;unique(hv.IngenuityCanonicalPathways)];
%    end
   FN.Pathways{z} = hv;
   clear hv
end

% create summary to plot pathway and GWAS overlap of CPA:
if (sum(contains(FN.Properties.VariableNames,'GSE'))>0)
    filesForOverlap = CreateFilesForOverlap(FN);
end


for z = 1: length(FN.Pathways)
   FN.Pathways{z} = FN.Pathways{z}(FN.Pathways{z}.logpvalue>-log10(0.05),:); % remove paths not enriched
   if z == 1 
       AllEnrichedPaths = unique(FN.Pathways{z}.IngenuityCanonicalPathways);
   else
       AllEnrichedPaths = [AllEnrichedPaths;unique(FN.Pathways{z}.IngenuityCanonicalPathways)];
   end
end


AllEnrichedPaths = table(AllEnrichedPaths,'variablenames',{'IngenuityCanonicalPathways'});
AllEnrichedPaths = varfun(@mean, AllEnrichedPaths,'groupingvariables',{'IngenuityCanonicalPathways'});
AllEnrichedPaths = sortrows(AllEnrichedPaths,'IngenuityCanonicalPathways','descend');

specialCharacters = {'Î±','Î²','Î³','Îº','Ïƒ'}; %{'α','β','γ','κ','σ'}
AllEnrichedPaths.IngenuityCanonicalPathways(contains(AllEnrichedPaths.IngenuityCanonicalPathways,specialCharacters));
AllEnrichedPaths = sortrows(AllEnrichedPaths,'GroupCount','descend');


%% get involved genes:
for z = 1 : length(FN.file)
    hv = FN.Pathways{z}(:,{'IngenuityCanonicalPathways','Molecules'});
    AllEnrichedPaths = outerjoin(AllEnrichedPaths,hv,'keys','IngenuityCanonicalPathways','mergekeys',1);
    clear hv
end
pos = find(contains(AllEnrichedPaths.Properties.VariableNames,'Molecules'));
AllEnrichedPaths.key = cellfun(@(x) sprintf('X%s',x),strtrim(cellstr(num2str([1:length(AllEnrichedPaths.GroupCount)]'))),'UniformOutput',0);
for z  = 1 : length(AllEnrichedPaths.GroupCount)
    AllEnrichedPaths.AllMolecules{z} = strjoin(table2array(AllEnrichedPaths(z,pos)),',');
    if contains(AllEnrichedPaths.AllMolecules{z},specialCharacters)
        error('gene names were compromised!')
    end
    AllEnrichedPaths.AllMolecules{z} = regexp(AllEnrichedPaths.AllMolecules{z},'\,','split');
    AllEnrichedPaths.AllMolecules{z}(strcmp(AllEnrichedPaths.AllMolecules{z},''))=[];
    AllEnrichedPaths.AllMolecules{z} = unique(AllEnrichedPaths.AllMolecules{z});
end
AllEnrichedPaths(:,pos)=[];
AllEnrichedPaths.pathsize = cellfun(@length, AllEnrichedPaths.AllMolecules);

%% get z-scores and p-values:
idx_joint = contains(cell_orgainzation_savename,'likeJoint');
idx_muscle = contains(cell_orgainzation_savename,'likeMuscle');

for z = 1 : length(AllEnrichedPaths.IngenuityCanonicalPathways)
    ZandP = NaN(length(cell_orgainzation),2);
    path = AllEnrichedPaths.IngenuityCanonicalPathways{z};
    for zz = 1 : length(cell_orgainzation)
        hv = FN.Pathways{strcmp(FN.key,cell_orgainzation{zz})};
        if sum(strcmp(hv.IngenuityCanonicalPathways,path))==0
            ZandP(zz,:) = [NaN,NaN];
        else
            ZandP(zz,1) = hv.logpvalue(strcmp(hv.IngenuityCanonicalPathways,path));
            ZandP(zz,2) = hv.zscore(strcmp(hv.IngenuityCanonicalPathways,path));
        end
    end
    ZandP(ZandP(:,2)==0,2)=NaN;
    if z == 1
        ZandPTotal = [table(repmat({path},size(ZandP,1),1),cell_orgainzation_savename,'variablenames',{'IngenuityCanonicalPathways','Cell'}),...
                       array2table(ZandP,'variablenames',{'logpval','zscore'}) ];
    else
        ZandPTotal = [ZandPTotal;[table(repmat({path},size(ZandP,1),1),cell_orgainzation_savename,'variablenames',{'IngenuityCanonicalPathways','Cell'}),...
                       array2table(ZandP,'variablenames',{'logpval','zscore'}) ]];
    end
    

    AllEnrichedPaths.M_activated_ratio(z) = sum(ZandP(idx_muscle,2)>0)/(2*sum(idx_muscle));
    AllEnrichedPaths.M_inhibited_ratio(z) = sum(ZandP(idx_muscle,2)<0)/(2*sum(idx_muscle));
    AllEnrichedPaths.M_noZ_ratio(z) = sum(isnan(ZandP(idx_muscle,2)) &    ZandP(idx_muscle,1)>-log10(0.05)    )/(2*sum(idx_muscle));
    AllEnrichedPaths.M_notSig_ratio(z) = sum(    isnan(ZandP(idx_muscle,1)) |    ZandP(idx_muscle,1)<-log10(0.05)    )/(2*sum(idx_muscle));
    
    AllEnrichedPaths.J_activated_ratio(z) = sum(ZandP(idx_joint,2)>0)/(2*sum(idx_joint));
    AllEnrichedPaths.J_inhibited_ratio(z) = sum(ZandP(idx_joint,2)<0)/(2*sum(idx_joint));
    AllEnrichedPaths.J_noZ_ratio(z) = sum(    isnan(ZandP(idx_joint,2)) &    ZandP(idx_joint,1)>-log10(0.05)    )/(2*sum(idx_joint));
    AllEnrichedPaths.J_notSig_ratio(z) = sum(    isnan(ZandP(idx_joint,1)) |    ZandP(idx_joint,1)<-log10(0.05)    )/(2*sum(idx_joint));    
    clear hv path ZandP
end 

% comparing proportions in cells/datasets where pathway is enriched.
AllEnrichedPaths.M_GlobalActivated = repmat(0,length(AllEnrichedPaths.key),1);
AllEnrichedPaths.M_GlobalInhibited = repmat(0,length(AllEnrichedPaths.key),1);
AllEnrichedPaths.M_GlobalNoDirection = repmat(0,length(AllEnrichedPaths.key),1);
AllEnrichedPaths.M_GlobalNotSignificant = repmat(0,length(AllEnrichedPaths.key),1);
idx = AllEnrichedPaths.M_activated_ratio>AllEnrichedPaths.M_inhibited_ratio & AllEnrichedPaths.M_activated_ratio>=AllEnrichedPaths.M_noZ_ratio & AllEnrichedPaths.M_activated_ratio>0;
AllEnrichedPaths.M_GlobalActivated(idx)=0.5;
clear idx
idx = AllEnrichedPaths.M_inhibited_ratio>AllEnrichedPaths.M_activated_ratio & AllEnrichedPaths.M_inhibited_ratio>=AllEnrichedPaths.M_noZ_ratio & AllEnrichedPaths.M_inhibited_ratio>0;
AllEnrichedPaths.M_GlobalInhibited(idx)=0.5;
clear idx
idx = AllEnrichedPaths.M_noZ_ratio>AllEnrichedPaths.M_activated_ratio & AllEnrichedPaths.M_noZ_ratio>AllEnrichedPaths.M_inhibited_ratio & AllEnrichedPaths.M_noZ_ratio>0;
AllEnrichedPaths.M_GlobalNoDirection(idx)=0.5;
clear idx
idx = AllEnrichedPaths.M_activated_ratio==AllEnrichedPaths.M_inhibited_ratio & AllEnrichedPaths.M_inhibited_ratio>=AllEnrichedPaths.M_noZ_ratio  & (AllEnrichedPaths.M_noZ_ratio>0 | AllEnrichedPaths.M_inhibited_ratio>0);
AllEnrichedPaths.M_GlobalNoDirection(idx)=0.5;
clear idx
idx = AllEnrichedPaths.M_activated_ratio==0 & AllEnrichedPaths.M_inhibited_ratio==0 & AllEnrichedPaths.M_noZ_ratio==0;
AllEnrichedPaths.M_GlobalNotSignificant(idx)=0.5;
clear idx

AllEnrichedPaths.J_GlobalActivated = repmat(0,length(AllEnrichedPaths.key),1);
AllEnrichedPaths.J_GlobalInhibited = repmat(0,length(AllEnrichedPaths.key),1);
AllEnrichedPaths.J_GlobalNoDirection = repmat(0,length(AllEnrichedPaths.key),1);
AllEnrichedPaths.J_GlobalNotSignificant = repmat(0,length(AllEnrichedPaths.key),1);

idx = AllEnrichedPaths.J_activated_ratio>AllEnrichedPaths.J_inhibited_ratio & AllEnrichedPaths.J_activated_ratio>=AllEnrichedPaths.J_noZ_ratio & AllEnrichedPaths.J_activated_ratio>0;
AllEnrichedPaths.J_GlobalActivated(idx)=0.5;
clear idx
idx = AllEnrichedPaths.J_inhibited_ratio>AllEnrichedPaths.J_activated_ratio & AllEnrichedPaths.J_inhibited_ratio>=AllEnrichedPaths.J_noZ_ratio & AllEnrichedPaths.J_inhibited_ratio>0;
AllEnrichedPaths.J_GlobalInhibited(idx)=0.5;
clear idx
idx = AllEnrichedPaths.J_noZ_ratio>AllEnrichedPaths.J_activated_ratio & AllEnrichedPaths.J_noZ_ratio>AllEnrichedPaths.J_inhibited_ratio & AllEnrichedPaths.J_noZ_ratio>0;
AllEnrichedPaths.J_GlobalNoDirection(idx)=0.5;
clear idx
idx = AllEnrichedPaths.J_activated_ratio==AllEnrichedPaths.J_inhibited_ratio & AllEnrichedPaths.J_inhibited_ratio>=AllEnrichedPaths.J_noZ_ratio & (AllEnrichedPaths.J_noZ_ratio>0 | AllEnrichedPaths.J_inhibited_ratio>0);
AllEnrichedPaths.J_GlobalNoDirection(idx)=0.5;
clear idx
idx = AllEnrichedPaths.J_activated_ratio==0 & AllEnrichedPaths.J_inhibited_ratio==0 & AllEnrichedPaths.J_noZ_ratio==0;
AllEnrichedPaths.J_GlobalNotSignificant(idx)=0.5;
clear idx

AllEnrichedPaths.activeTissue = repmat({''},length(AllEnrichedPaths.AllMolecules),1);
AllEnrichedPaths.activeTissue(AllEnrichedPaths.J_GlobalActivated==0.5)={'activated'};
AllEnrichedPaths.activeTissue(AllEnrichedPaths.J_GlobalInhibited==0.5)={'inhibited'};
AllEnrichedPaths.activeTissue(AllEnrichedPaths.J_GlobalNoDirection==0.5)={'no_direction'};
AllEnrichedPaths.activeTissue(AllEnrichedPaths.J_GlobalNotSignificant==0.5)={'not_sig'};

AllEnrichedPaths.inactiveTissue = repmat({''},length(AllEnrichedPaths.AllMolecules),1);
AllEnrichedPaths.inactiveTissue(AllEnrichedPaths.M_GlobalActivated==0.5)={'activated'};
AllEnrichedPaths.inactiveTissue(AllEnrichedPaths.M_GlobalInhibited==0.5)={'inhibited'};
AllEnrichedPaths.inactiveTissue(AllEnrichedPaths.M_GlobalNoDirection==0.5)={'no_direction'};
AllEnrichedPaths.inactiveTissue(AllEnrichedPaths.M_GlobalNotSignificant==0.5)={'not_sig'};
