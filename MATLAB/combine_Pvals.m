function res = combine_Pvals(FishertestR)
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
    res = CombinedPval_Disease;

end