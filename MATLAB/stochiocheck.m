lida_petrinet = sbioloadproject('LIDA.sbproj').m1;
lida_mine = buildmodel('lida','explore');

S1 = full(getstoichmatrix(lida_petrinet));
S2 = full(getstoichmatrix(lida_mine));

swap = [];
for i = 1:length(S1)
for j = 1:length(S1)
if S1(j,:) == S2(i,:)
swap = [swap;[i j]]; %#ok<AGROW>
end
end
end

S1 = S1(swap(:,2),:);

if S1 == S2
sprintf('ok')
else
sprintf('not ok')
end


ODE1 = getequations(lida_petrinet)
ODE2 = getequations(lida_mine)