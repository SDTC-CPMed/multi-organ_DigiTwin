function [res1 , res2] = preprocess(gene_info, AllDatasets)

    
    % translate all genes into gene symbols
    for z = 1 : length(AllDatasets.deg_file_name)
        hv = readtable(join(['../data/AllDEGfilesMovedToOneFolder/', AllDatasets.deg_file_name{z}]));
        if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE81071_DLE_vs_control.csv','GSE81071_SCLE_vs_control.csv','GSE95065_SSC_skin.csv'}))==1
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ENTREZ_GENE_ID'))={'Gene_ID'};
        end
        if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE148810_juvenile myositis_skin_1.csv','GSE148810_Nonlesional skin 6 vs HC 8_JM.csv'}))==1
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ORF'))={'Gene_symbol'};
        end
        if strcmp(AllDatasets.deg_file_name{z},'GSE32591_glomer_vs_contol_LN.csv')
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'Gene_Symbol'))={'Gene_symbol'};
        end
        if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE181318_skin_psoriatic 3 vs control3.csv','GSE66413_Pancreatic lymph nodes 13_ T1D vs healthy 3.csv'}))
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'GENE_SYMBOL'))={'Gene_symbol'};
        end
        if strcmp(AllDatasets.deg_file_name{z},'GSE176510_Sj”gren syndrome_keratocojunctivitis scicca.csv')
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ID'))={'Gene_symbol'};
        end
        if strcmp(AllDatasets.deg_file_name{z},'GSE81292_SSc_ILD_lung.csv')
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ORF'))={'Gene_ID'};
        end
        if strcmp(AllDatasets.deg_file_name{z},'GSE148810_childhood_onset_lupus_cSLE_skin.csv')
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'ORF'))={'Gene_symbol'};
            hv.Properties.VariableNames(strcmp(hv.Properties.VariableNames,'SPOT_ID'))={'Gene_ID'};
        end 

        if sum(strcmp(AllDatasets.deg_file_name{z},{'GSE75214_inactive_vs_normal_CD_16_11_ileum.csv','GSE148810_juvenile myositis_skin_1.csv',...
                'GSE148810_Nonlesional skin 6 vs HC 8_JM.csv','GSE32591_glomer_vs_contol_LN.csv','GSE181318_skin_psoriatic 3 vs control3.csv',...
                'GSE176510_Sj”gren syndrome_keratocojunctivitis scicca.csv','GSE75214_inactive_vs_normal_23_11_UC_colon.csv','GSE66413_Pancreatic lymph nodes 13_ T1D vs healthy 3.csv'}))==1
        % For cases of genes with alternative names

            % collect all measured genes
            genesfound = unique(hv.Gene_symbol(~strcmp(hv.Gene_symbol,'')));
            for i = 1:length(genesfound)
                genes = regexp(genesfound(i),'///','split'); 
                for j = 1:length(genes{1,1})
                    if ismember(genes{1,1}(j), gene_info.Symbol) 
                        genesfound(i) = genes{1,1}(j);
                    end
                end
            end
            genesfound(~ismember(genesfound,gene_info.Symbol),:)=[]; 
            AllDatasets.measured_genes{z} = genesfound;

            % collect all differentially expressed genes
            hv(hv.adj_P_Val>=0.05,:)=[];
            genesfound = unique(hv.Gene_symbol(~strcmp(hv.Gene_symbol,'')));
            for i = 1:length(genesfound)
                genes = regexp(genesfound(i),'///','split');
                for j = 1:length(genes{1,1})
                    if ismember(genes{1,1}(j), gene_info.Symbol)
                        genesfound(i) = genes{1,1}(j);
                    end
                end
            end
            genesfound(~ismember(genesfound,gene_info.Symbol),:)=[];
            AllDatasets.DEG_genes{z}  = unique(genesfound);

        else 
            % collect all measured genes
            genesfound = unique(hv.Gene_ID(~isnan(hv.Gene_ID)));
            genesfound(~ismember(genesfound,gene_info.GeneID),:)=[];
            AllDatasets.measured_genes{z}  = unique(gene_info.Symbol(ismember(gene_info.GeneID,genesfound)));

            % collect all differentially expressed genes        
            hv(hv.adj_P_Val>=0.05,:)=[];
            genesfound(~ismember(genesfound,hv.Gene_ID),:)=[];
            AllDatasets.DEG_genes{z}  = unique(gene_info.Symbol(ismember(gene_info.GeneID,genesfound)));
        end
        clear genesfound    
    end

    %% Restructuring the dataset

    AllDatasets.broad_label = repmat({''},length(AllDatasets.key),1);
    AllDatasets.broad_label(strcmp(AllDatasets.active_USE,'yes')) = cellfun(@(x) sprintf('%s_active',x),AllDatasets.disease_short(strcmp(AllDatasets.active_USE,'yes')),'UniformOutput',0);
    AllDatasets.broad_label(strcmp(AllDatasets.active_USE,'no')) = cellfun(@(x) sprintf('%s_inactive',x),AllDatasets.disease_short(strcmp(AllDatasets.active_USE,'no')),'UniformOutput',0);

    AllDatasets.DEG=[];
    AllDatasets.group = repmat({'Inflamed'},length(AllDatasets.DEG_genes),1);
    AllDatasets.group(strcmp(AllDatasets.active_USE,'no'))={'Noninflamed'};
    AllDatasets.Properties.VariableNames(strcmp(AllDatasets.Properties.VariableNames,'DEG_genes'))={'DEG'};
    AllDatasets.Properties.VariableNames(strcmp(AllDatasets.Properties.VariableNames,'disease_short'))={'cell_dataset'};
    AllDatasets.key = cellfun(@(x,y) sprintf('%s***%s',x,y),AllDatasets.file_name,AllDatasets.group,'UniformOutput',0);


    broad_label = AllDatasets.broad_label;
    broad_label(10) = {'lupus_active'};
    broad_label(13) = {'lupus_active'};
    broad_label(18) = {'lupus_active'};
    broad_label(30) = {'lupus_active'};
    AllDatasets.broad_label = broad_label;

    Fn = AllDatasets(:,{'key','file_name','group','DEG', 'broad_label'});
    res1 = Fn;
    res2 = AllDatasets;
end