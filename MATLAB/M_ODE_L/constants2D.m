NA = 6.023e23;   
D = 9.32e-10;   %filippov2004 taken in dm^2/s
gamma = 0.57722;
R = 2e-8;       %recradius 2nm 
SF = 1/200;

phi = pi*R^2*1e-4;
logphi = log(phi);
den = log(4*D*2e6/(R^2)) - 2*gamma;

%% RateConsts

kf_or1_tleft    = (4*pi*D*NA*SF/den);        %2D on-rate for oligomer1+
kf_or2_tleft    = 2e7;                       % on-rate for oligomer2+
kf_ol1_tright   = (4*pi*D*NA*SF/den);        %2D on-rate for oligomer1-
kf_ol2_tright   = 2e7;                       % on-rate for oligomer2-
                             
kb_or1_tleft   = 7.55;     % off-rate for oligomer1+
kb_or2_tleft   = 6.7950e4; % off-rate for oligomer2+
kb_ol1_tright  = 15.824;   % off-rate for oligomer1-
kb_ol2_tright  = 7.156e4;  % off-rate for oligomer2- 

kf_repair  = 2e-2;                 % reparation rate from nicked-duplex to double strand in the left one
kf_duplex  = (4*pi*D*NA*SF/den);   %2D on-rate for templates
kb_duplex = 1.3776e-4;             % off-rate for templates

% Background reaction constants
kf_orl_1  = (4*pi*D*NA*SF/den);    %2D on-rate for destabilizing strands 1
kf_orl_2  = 2e7;                   % on-rate for destabilizing strands 2
kb_orl_1 = 2e1;                  % off-rate for destabilizing strands 1
kb_orl_2 = 1e5;                    % off-rate for destabilizing strands 2
kpara = 0;                              % Ligation of the two O1+- and O2+- 

