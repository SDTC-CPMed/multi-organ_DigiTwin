function [res1, res2] = count_significant(FishertestR)
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
    res1 = count_Inf;
    res2 = count_Noninf;

end