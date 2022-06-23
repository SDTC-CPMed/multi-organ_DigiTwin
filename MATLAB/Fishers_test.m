function res = Fishers_test(SP,UR,Fn)

    FPval = cell(length(SP.SP),1);
    FOR = cell(length(SP.SP),1);


    uUR = unique(UR.UR);
    for p = 1 : length(SP.SP)
        for cd = 1 : length(Fn.key)               
                spgenes = SP.AllMolecules{p}(ismember(SP.AllMolecules{p},Fn.DEG{cd}));
                URhv = UR(strcmp(UR.key,Fn.key{cd}),:);

                for ur = 1 : length(uUR)
                    if sum(strcmp(URhv.UR,uUR{ur}))
                        if length(unique(Fn.DEG{cd})) ~= length(Fn.DEG{cd})
                            error = p;
                        end

                        urds = URhv.DS{strcmp(URhv.UR,uUR{ur})}; 

                        a = sum(ismember(spgenes,urds));
                        b = sum(~ismember(spgenes,urds));
                        c = sum(~ismember(urds,spgenes));
                        d = sum(~ismember(Fn.DEG{cd},[spgenes,urds]));

                        [~, FPval{p}(cd,ur)] = fishertest([a b; c d],'tail','right');
                        FOR{p}(cd,ur) = (a*d)/(b*c);
                    else
                        FPval{p}(cd,ur) = 1;
                        FOR{p}(cd,ur) = 0;
                    end

                    clear urds a b c d
                end
                clear spgenes
        end
    end

    % Summarize results
    FishertestR.cellsdes = SP;
    FishertestR.rowdes = Fn;
    FishertestR.coldes_URs = uUR;
    FishertestR.moreURinfo = UR;
    FishertestR.Pval = FPval;
    FishertestR.OR = FOR;
    
    res = FishertestR;
end