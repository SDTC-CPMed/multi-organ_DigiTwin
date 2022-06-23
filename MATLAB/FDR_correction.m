function res = FDR_correction(CombinedOverInfandNoninf)
    for i = 1:length(CombinedOverInfandNoninf.CombP_Inf(:,1))
        FDR_Inf(i,:) = mafdr(CombinedOverInfandNoninf.CombP_Inf(i,:), 'BHFDR',true);
        FDR_Noninf(i,:) = mafdr(CombinedOverInfandNoninf.CombP_Noninf(i,:), 'BHFDR',true);
        FDR_All(i,:) = mafdr(CombinedOverInfandNoninf.CombP_all(i,:), 'BHFDR',true);

    end

    CombinedOverInfandNoninf.FDR_Inf = FDR_Inf;
    CombinedOverInfandNoninf.FDR_Noninf = FDR_Noninf;
    CombinedOverInfandNoninf.FDR_All = FDR_All;

    res = CombinedOverInfandNoninf;
end