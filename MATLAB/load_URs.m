function res = load_URs(AllDatasets)
    AllDatasets = AllDatasets(:,{'file_name','key','cell_dataset','tissue','gse','active_USE','UR', 'measured_genes','DEG','group'});

    for z = 1 : length(AllDatasets.key)
        AllDatasets.UR{z}(AllDatasets.UR{z}.p_valueofoverlap>0.05,:)=[]; 
        hv = [table(repmat(AllDatasets.file_name(z),length(AllDatasets.UR{z}.upstreamregulator),1),...
            repmat(AllDatasets.group(z),length(AllDatasets.UR{z}.upstreamregulator),1),'variablenames',{'file_name','group'}),...
            AllDatasets.UR{z}(:,{'upstreamregulator','targetmoleculesindataset'})];
        if z == 1
            allURhv = hv;
        else        
            allURhv = [allURhv;hv];
        end
    end
    allURhv.targetmoleculesindataset =  regexp(allURhv.targetmoleculesindataset,'\,','split');
    % Check if DSs are among DEGs
    for i = 1:length(allURhv.targetmoleculesindataset)
        DEGs = AllDatasets.DEG{strcmp(AllDatasets.file_name, allURhv.file_name(i))};
        allURhv.targetmoleculesindataset{i,1}(~ismember(allURhv.targetmoleculesindataset{i,1},DEGs)) = [];
    end

    allURhv.Properties.VariableNames(strcmp(allURhv.Properties.VariableNames,'targetmoleculesindataset'))={'DS'};
    allURhv.Properties.VariableNames(strcmp(allURhv.Properties.VariableNames,'upstreamregulator'))={'UR'};
    allURhv.key  = cellfun(@(x,y) sprintf('%s***%s',x,y),allURhv.file_name,allURhv.group,'UniformOutput',0);
    UR = allURhv(:,{'UR','key','file_name','group','DS'});
    res = UR;
end