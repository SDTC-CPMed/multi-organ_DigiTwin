function res = format_output(CombinedOverInfandNoninf,CombinedPval_Disease)
    CombinedPval  = table(repmat(CombinedOverInfandNoninf.rows_SPs,size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    repelem(CombinedOverInfandNoninf.cols_URs,size(CombinedOverInfandNoninf.CombP_Inf,1)),...
                    reshape(CombinedOverInfandNoninf.CombP_Inf,size(CombinedOverInfandNoninf.CombP_Inf,1)*size(CombinedOverInfandNoninf.CombP_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_Noninf,size(CombinedOverInfandNoninf.CombP_Noninf,1)*size(CombinedOverInfandNoninf.CombP_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.CombP_all,size(CombinedOverInfandNoninf.CombP_all,1)*size(CombinedOverInfandNoninf.CombP_all,2),1),...
                    reshape(CombinedOverInfandNoninf.FDR_Inf,size(CombinedOverInfandNoninf.FDR_Inf,1)*size(CombinedOverInfandNoninf.FDR_Inf,2),1),...
                    reshape(CombinedOverInfandNoninf.FDR_Noninf,size(CombinedOverInfandNoninf.FDR_Noninf,1)*size(CombinedOverInfandNoninf.FDR_Noninf,2),1),...
                    reshape(CombinedOverInfandNoninf.FDR_All,size(CombinedOverInfandNoninf.FDR_All,1)*size(CombinedOverInfandNoninf.FDR_All,2),1),...
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
                    'variablenames',{'SP','UR','CombinedP_Inflamed','CombinedP_Noninflamed', 'CombinedP_All', 'FDR_Inflamed', 'FDR_Noninflamed', 'FDR_All', 'count_Inflamed', 'count_Noninflamed','AD_active','AD_inactive','CD_active','CD_inactive','JM_active','JM_inactive','PSO_active','PSO_inactive','RA_active','SS_active','SSc_active','UC_active','UC_inactive','at_risk_T1D_inactive','lupus_active'});

    res = CombinedPval;
end