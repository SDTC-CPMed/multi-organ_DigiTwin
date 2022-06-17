function Jaccard = CalculateJaccardIndex(AllEnrichedPaths)

Jaccard = NaN(length(AllEnrichedPaths.AllMolecules));
for z = 1 : length(AllEnrichedPaths.AllMolecules)
   for zz = 1 :  length(AllEnrichedPaths.AllMolecules)
       Jaccard(z,zz) = sum(ismember(AllEnrichedPaths.AllMolecules{z},AllEnrichedPaths.AllMolecules{zz}))/length(unique([AllEnrichedPaths.AllMolecules{z},AllEnrichedPaths.AllMolecules{zz}]));
   end
end
% issymmetric(Jaccard)
Jaccard = array2table(Jaccard,'variablenames',AllEnrichedPaths.key,'rownames',AllEnrichedPaths.key);
