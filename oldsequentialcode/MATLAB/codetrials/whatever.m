clear all
%#ok<*NOPTS>

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

sbioloadproject LIDA_noback_versionuniquecompartment.sbproj


[Sigma,Species,Reactions] = getstoichmatrix(m1);
Sigma = full(Sigma);

swap = [];

for i = 1:length(Sigma)
for j = 1:length(Sigma)
if Sigma(j,:) == Stochio(i,:)
swap = [swap;[i j]]; %#ok<AGROW>
end
end
end


%La matrice di simbiology permutando le righe in modo da matchare il mio
%ordine mentale delle reazioni

csobj = getconfigset(m1);
set(csobj,'Stoptime',4e4);
set(csobj,'TimeUnits','hour');
set(csobj.SolverOptions,'AbsoluteTolerance',1e-16);
get(csobj.SolverOptions,'AbsoluteTolerance');
set(csobj.SolverOptions,'RelativeTolerance',5e-13);

[t,x,n] = sbiosimulate(m1);

% x = x(swap(:,2),:);
% n = n(swap(:,2),:);

c = x(:,swap(:,2));

plot(t,c)
legend("Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L")

names = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];
[mxx,ixx] = max(c);

step = 4e4/4496
check = 6;

tmax = t(ixx(check))

figure(check)
hold on
plot(t,c(:,check))
plot(tmax,mxx(check),'r+')


