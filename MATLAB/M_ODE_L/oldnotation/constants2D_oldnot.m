% 
% NA = 6.023e23;   
% D = 9.32e-10;   %filippov2004 in dm^2/s
% gamma = 0.57722;
% R = 2e-8;       %recradius 2nm 
% SF = 1/200;
% phi = pi*R^2*1e-4;
% den = log(4*D*2000/(R^2)) - 2*gamma;

% ===================================================  2D Constants
%% set2

kpluso1tplus    = paramvalue(1);             %2D on-rate for oligomer1+
kpluso2tplus    = paramvalue(3);                       % on-rate for oligomer2+
kpluso1tminus   = paramvalue(1);             %2D on-rate for oligomer1-
kpluso2tminus   = paramvalue(3);                       % on-rate for oligomer2-
                             
kminuso1tplus   = paramvalue(4);                    % off-rate for oligomer1+
kminuso2tplus   = paramvalue(5);                       % off-rate for oligomer2+
kminuso1tminus  = paramvalue(4);                     % off-rate for oligomer1-
kminuso2tminus  = paramvalue(5);                        % off-rate for oligomer2-

klplus  = paramvalue(7);                              % ligation rate from nicked-duplex to double strand in the plus one
klminus = paramvalue(7);                              % ligation rate from nicked-duplex to double strand in the minus one
ktplus  = paramvalue(1);                        %2D on-rate for templates
ktminus = paramvalue(6);                             % off-rate for templates

% Background reaction constants
kpluso1tpm  = paramvalue(1);                 %2D on-rate for destabilizing strands 1
kpluso2tpm  = paramvalue(2);                           % on-rate for destabilizing strands 2
kminuso1tpm = paramvalue(4);                        % off-rate for destabilizing strands 1
kminuso2tpm = paramvalue(5);                            % off-rate for destabilizing strands 2

 

%% set1
kpluso1tplus    = kf_2d.value;             %2D on-rate for oligomer1+
kpluso2tplus    = kf_mix.value;                       % on-rate for oligomer2+
kpluso1tminus   = kf_2d.value;             %2D on-rate for oligomer1-
kpluso2tminus   = kf_mix.value;                       % on-rate for oligomer2-
                             
kminuso1tplus   = kb_oligo_surf.value;                    % off-rate for oligomer1+
kminuso2tplus   = kb_oligo_bulk.value;                       % off-rate for oligomer2+
kminuso1tminus  = kb_oligo_surf.value;                     % off-rate for oligomer1-
kminuso2tminus  = kb_oligo_bulk.value;                        % off-rate for oligomer2-

klplus  = 2e-2;                              % ligation rate from nicked-duplex to double strand in the plus one
klminus = 2e-2;                              % ligation rate from nicked-duplex to double strand in the minus one
ktplus  = kf_2d.value;                        %2D on-rate for templates
ktminus = kb_duplex.value;                             % off-rate for templates

% Background reaction constants
kpluso1tpm  = kf_2d.value;                 %2D on-rate for destabilizing strands 1
kpluso2tpm  = kf_3d.value;                           % on-rate for destabilizing strands 2
kminuso1tpm = kb_oligo_surf.value;                        % off-rate for destabilizing strands 1
kminuso2tpm = kb_oligo_bulk.value;                            % off-rate for destabilizing strands 2

 

% kpluso1tplus    = 4*pi*D*NA*SF/den;             %2D on-rate for oligomer1+
% kpluso2tplus    = 2e7;                       % on-rate for oligomer2+
% kpluso1tminus   = 4*pi*D*NA*SF/den;             %2D on-rate for oligomer1-
% kpluso2tminus   = 2e7;                       % on-rate for oligomer2-
%                              
% kminuso1tplus   = 7.16e5;                    % off-rate for oligomer1+
% kminuso2tplus   = 158.24;                       % off-rate for oligomer2+
% kminuso1tminus  = 6.8e5;                     % off-rate for oligomer1-
% kminuso2tminus  = 75.46;                        % off-rate for oligomer2-
% 
% klplus  = 2e-2;                              % ligation rate from nicked-duplex to double strand in the plus one
% klminus = 2e-2;                              % ligation rate from nicked-duplex to double strand in the minus one
% ktplus  = 4*pi*D*NA*SF/den;                        %2D on-rate for templates
% ktminus = 0.3767;                             % off-rate for templates
% 
% % Background reaction constants
% kpluso1tpm  = 4*pi*D*NA*SF/den;                 %2D on-rate for destabilizing strands 1
% kpluso2tpm  = 2e7;                           % on-rate for destabilizing strands 2
% kminuso1tpm = 5.62e4;                        % off-rate for destabilizing strands 1
% kminuso2tpm = 163.8;                            % off-rate for destabilizing strands 2


