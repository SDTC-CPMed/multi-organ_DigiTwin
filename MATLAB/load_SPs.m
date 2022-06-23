function res = load_SPs(path_cpa_program1, path_cpa_program2)

    P1 = readtable(path_cpa_program1);
    P2 = readtable(path_cpa_program2);

    P1.subprograms = cellfun(@(x) sprintf('1.%s',x), strtrim(cellstr(num2str(P1.subclusters))),'UniformOutput',0);
    P2.subprograms = cellfun(@(x) sprintf('2.%s',x), strtrim(cellstr(num2str(P2.subclusters))),'UniformOutput',0);
    P = [P1(:,{'clustersGlobal','subprograms','AllMolecules'});P2(:,{'clustersGlobal','subprograms','AllMolecules'})];
    Programs = table(unique(P.subprograms),'variablenames',{'Program'});

    for z = 1 : length(Programs.Program)
        Programs.AllMolecules{z} = unique(regexp(strjoin(P.AllMolecules(strcmp(P.subprograms,Programs.Program{z})),','),'\,','split'));
    end
    Programs.Program(length(Programs.Program)+1)={'P1'};
    Programs.AllMolecules{end} = unique(regexp(strjoin(P.AllMolecules(P.clustersGlobal==1),','),'\,','split'));

    Programs.Program(length(Programs.Program)+1)={'P2'};
    Programs.AllMolecules{end} = unique(regexp(strjoin(P.AllMolecules(P.clustersGlobal==2),','),'\,','split'));


    SP=Programs;
    SP.Properties.VariableNames(strcmp(SP.Properties.VariableNames,'Program'))={'SP'};
    res = SP;
end