
syms E
syms S
syms C
syms P
syms kf1
syms km1
syms kf2
syms edot
syms sdot
syms cdot
syms pdot


Species = [E S C P];
Sigma = [-1 1 1;-1 1 0;1 -1 -1;0 0 1];
M = [E*S,C,C];
K = [kf1,km1,kf2];
Gamma = diag(M)*diag(K);
O = Sigma*Gamma;

OD = sum(O,2);

dots = [-edot; -sdot; -cdot; -pdot];

O = [O dots];

null(O)
orth(O)

% s0 = 10;
% e0 = 1;
% 
% dq = @(t,y)[-10*y(1)*(e0-y(2))+8*y(2);10*y(1)*(e0-y(2))-8*y(2)-12*y(2);12*y(2)];
% [t,sol] = ode45(dq,(0:1e-4:2),[10,0,0]);
% 
% plot(t,sol)
% plot(t,e0-sol(:,2))
