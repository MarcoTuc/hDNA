kf_2d = 1.5e8;

kpluso1tplus    = kf_2d;             %2D on-rate for oligomer1+
kpluso2tplus    = 2e7;                       % on-rate for oligomer2+
kpluso1tminus   = kf_2d;             %2D on-rate for oligomer1-
kpluso2tminus   = 2e7;                       % on-rate for oligomer2-
                             
kminuso1tplus   = 71600;                    % off-rate for oligomer1+
kminuso2tplus   = 158;                       % off-rate for oligomer2+
kminuso1tminus  = 68000;                     % off-rate for oligomer1-
kminuso2tminus  = 75.5;                        % off-rate for oligomer2-

klplus  = 2e-2;                              % ligation rate from nicked-duplex to double strand in the plus one
klminus = 2e-2;                              % ligation rate from nicked-duplex to double strand in the minus one
ktplus  = kf_2d;                        %2D on-rate for templates
ktminus = 0.37;                             % off-rate for templates

% Background reaction constants
kpluso1tpm  = kf_2d;                 %2D on-rate for destabilizing strands 1
kpluso2tpm  = 2e7;                           % on-rate for destabilizing strands 2
kminuso1tpm = 56000;                        % off-rate for destabilizing strands 1
kminuso2tpm = 163.8;                            % off-rate for destabilizing strands 2

model2D_oldnot