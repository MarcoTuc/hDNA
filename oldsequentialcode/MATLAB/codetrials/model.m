clear all


sbioloadproject LIDA_noback.sbproj


[M,objSpecies] = getstoichmatrix(m1);

S = rref(M);


A = [   -1 0;
        -1 0;
        1 -1;
        0 -1;
        0 1]

R = [   -1 -1 1 0 0;
        -1 -1 1 0 0;
        1 1 -2 -1 0;
        0 0 -1 -1 0;
        0 0 1 1 0;]

U = pinv(A)*R


S1 = [-1 0;-1 0;1 -1;0 -1;0 1];
U1 = [5 5 -4 0 0;0 0 3 3 0];

syms e real
syms s real
syms c real 
syms a real
syms p real
syms k1 real
syms km1 real
syms k2 real
syms X
syms Y
syms B
syms C
syms a
syms km2

% Us = [kf*e*s kf*e*s -kb*c 0 0;0 0 k2*c*a k2*c*a 0];


A = [   1, 0;
        1, 0;
        -1, 1;
        0, 1;
        0, -1]



% transpose(A)*A
% A*transpose(A)

clear S G

S = [-2, 0;1, -1;0, 1;0, -1;0, 1]
G = [2*a*k1, -km1, 0, 0, 0; 0, k2*X, -Y*km2, B*k2, -km2*C]

J = S*G


Q = [k1 0 0; 0 km1 0; 0 0 k2];
E = [e*s  0 0;0 c 0; 0 0 c];
F = Q*E




 