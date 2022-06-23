function res = Fishers_method(FishertestR)
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


    % Create a summary file where we store important results
    CombinedOverInfandNoninf.cols_URs = FishertestR.coldes_URs;
    CombinedOverInfandNoninf.rows_SPs = FishertestR.cellsdes.SP;
    CombinedOverInfandNoninf.CombP_Inf = CombinedPval_Inflamed;
    CombinedOverInfandNoninf.CombP_Noninf = CombinedPval_Noninflamed;
    CombinedOverInfandNoninf.CombP_all = CombinedPval_all;
    
    res = CombinedOverInfandNoninf;
end