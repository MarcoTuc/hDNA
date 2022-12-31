clear all


sbioloadproject LIDA_noback_versionuniquecompartment.sbproj


[M,Species,Reactions] = getstoichmatrix(m1);

S = rref(M);
M = full(M)
Species;
Reactions;
sprank(M)
sprank(S)

rowNames = Species
colNames = Reactions
sTable = array2table(M,'RowNames',rowNames,'VariableNames',colNames)

idx = [1 2 3 5 4 6 7 8 9 10 11 12 13 14 15];

% S = M(idx,:)
% 
% rowNames = Species
% colNames = 
% sTable = array2table(S,'RowNames',rowNames,'VariableNames',colNames)

[W, pivots] = rref(M);
rank(full(M))
rank(rref(M))

[L U P] = lu(M)
rank(rref(L))
rank(U)
rref(L)
rref(U)
P
