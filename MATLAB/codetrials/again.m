clear all
%#ok<*NOPTS>

sbioloadproject LIDA_noback_versionuniquecompartment.sbproj


[Sigma,Species,Reactions] = getstoichmatrix(m1);

%Basic routine for getting linearly combined reactions
Sigma = full(Sigma);
Sig = Sigma';
[reducedSigma, Pivots] = rref(Sigma);

basis = transpose(Sig(Pivots,:));

reaction_4 = Sig(4,:);
x_4 = basis\reaction_4';

%% Trying it out on the ODE matrix

Stochio = [
%     1  2  3  4  5  6  7  8  9 10 11 12 13 Tl
    [-1,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 R1
    [-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0,-1, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 R2
    [ 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 TR1
    [ 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 TR2
    [ 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 Nl
    [ 0, 0, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 Tr
    [ 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 L1
    [ 0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 0,-1, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 L2
    [ 0, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0, 0,-1]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 TL1
    [ 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 TL2
    [ 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 Nr
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 D
    [ 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 0, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 O1
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
%     1  2  3  4  5  6  7  8  9 10 11 12 13 O2
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
]; 

rowNames = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];
colNames = Reactions;
sTable = array2table(Stochio,'RowNames',rowNames,'VariableNames',colNames);
[IR, Pivots] = rref(Stochio);
rank(Stochio);

rowNames = Species;
colNames = Reactions;
sTable = array2table(Sigma,'RowNames',rowNames,'VariableNames',colNames);
swap = [];

for i = 1:length(Sigma)
for j = 1:length(Sigma)
if Sigma(j,:) == Stochio(i,:)
swap = [swap;[i j]]; %#ok<AGROW>
end
end
end
swap

%La matrice di simbiology permutando le righe in modo da matchare il mio
%ordine mentale delle reazioni
swappata = Sigma(swap(:,2),:)

%% ora cerco di essere piu ordinato
S = swappata

syms Tl R1 R2 TR1 TR2 Nl Tr L1 L2 TL1 TL2 Nr D O1 O2 
syms kf_TR1 kb_TR1 kf_TR2 kb_TR2 kf_TL1 kb_TL1 kf_TL2 kb_TL2 kf_R kf_D kb_D kf_O1 kb_O1 kf_O2 kb_O2

% K = diag([kf_TR1 kb_TR1 kf_T])
% T = diag([Tl*R1 TR1 Tl*R2 TR2 TR1*R2 Nl TR2*R1 Nl Nl D D Tl*Tr Tr*L1 TL1 Tr*L2 TL2 TL1*L2 Nr TL2*L1 Nr Nr D R1*L1 O1 R1*L2 O2])

reactants = transpose([Tl,R1; Tl,R2; TR1,R2; TR2,R1; Nl,1; D,1; Tr,L1; Tr,L2; TL1,L2; TL2,L1; Nr,1; R1,L1; R2,L2])
products = transpose([TR1,TR2]);

Sf = prod(reactants(:,:))
Kf = [kf_TR1,kf_TR2,kf_TR2,kf_TR1,kf_R,kb_D,kf_TL1,kf_TL2,kf_TL2,kf_TL1,kf_R,kf_O1,kf_O2]

Sb = [TR1,TR2,Nl,Nl,D,Tl*Tr,TL1,TL2,Nr,Nr,D,O1,O2]
Kb = [kb_TR1,kb_TR2,kb_TR2,kb_TR1,0,kb_D,kb_TL1,kb_TL2,kb_TL2,kb_TL1,0,kb_O1,kb_O2]
GammaF = diag(Sf)*diag(Kf)
GammaB = diag(Sb)*diag(Kb)

Forward = S*GammaF
Backward = -S*GammaB

FOM = [Forward Backward]

% R = [T K]

[H] = rref(FOM)

