
%% ===================================================  Physically realistic

if kset == "realistic"

        kf_TR1.value   = 2e7;    %2D Template_Left + Oligomer_Right_1
        kf_TR1.name    = 'kf_TR1';
        kf_TR2.value   = 2e7;     %3D Template_Left + Oligomer_Right_2
        kf_TR2.name    = 'kf_TR2';
        kf_TL1.value   = 2e7;    %2D Template_Right + Oligomer_Left_1
        kf_TL1.name    = 'kf_TL1';
        kf_TL2.value   = 2e7;     %3D Template_Right + Oligomer_Left_2
        kf_TL2.name    = 'kf_TL2';

        kb_TR1.value   = 7.16e5;  %2D Template_Left + Oligomer_Right_1
        kb_TR1.name    = 'kb_TR1';
        kb_TR2.value   = 158.24;  %3D Template_Left + Oligomer_Right_2
        kb_TR2.name    = 'kb_TR2';
        kb_TL1.value   = 6.8e5;    %2D Template_Right + Oligomer_Left_1
        kb_TL1.name    = 'kb_TL1';
        kb_TL2.value   = 75.46;    %3D Template_Right + Oligomer_Left_2
        kb_TL2.name    = 'kb_TL2';

        kf_lig_left.value  = 2e-2;    % Nicked_Left  -> Duplex
        kf_lig_left.name   = 'kf_lig_left';

        kf_lig_right.value = 2e-2;    % Nicked_Right -> Duplex
        kf_lig_right.name  = 'kf_lig_right';
        kf_duplex.value    = 2e7;    % Simplexes -> Duplex
        kf_duplex.name     = 'kf_duplex';
        kb_duplex.value    = 0.3767;  % Duplex -> Simplexes
        kb_duplex.name     = 'kb_duplex';

        % Background reaction constants
        kf_O1.value  = 2e7;      %2D R1 + L1 -> O1
        kf_O1.name   = 'kf_O1';
        kf_O2.value  = 2e7;       %3D R2 + L2 -> O2 
        kf_O2.name   = 'kf_O2';
        kb_O1.value  = 5.62e4;     %2D O1 -> R1 + L1
        kb_O1.name   = 'kb_O1';
        kb_O2.value  = 163.8;      %3D O2 -> R2 + L2 
        kb_O2.name   = 'kb_O2';

        % not happening in a twodimensional LIDA
        % kpara = 2e-2;                                 % Ligation of the two O1+- and O2+- 

        %dummy rates for reversible
        kb_lig_left.value  = 0;
        kb_lig_left.name   = 'kb_lig_left';
        kb_lig_right.value = 0;
        kb_lig_right.name  = 'kb_lig_right';

parameters = [kf_TR1,kf_TR2,kf_TL1,kf_TL2,kb_TR1,kb_TR2,kb_TL1,kb_TL2,kf_lig,...
              kf_duplex,kb_duplex,kf_O1,kf_O2,kb_O1,kb_O2,kb_lig];

end

%% ===================================================  stochastic rate constants
if kset == "explore-stoch"

Na = 6.023e23;

k_lig = 2e-2 / Na;
kf_2d = 1.5e8 / Na;
kf_3d = 2e7 / Na;

        kf_TR1.value   = kf_2d;    %2D Template_Left + Oligomer_Right_1
        kf_TR1.name    = 'kf_TR1';
        kf_TR2.value   = kf_3d;     %3D Template_Left + Oligomer_Right_2
        kf_TR2.name    = 'kf_TR2';
        kf_TL1.value   = kf_2d;    %2D Template_Right + Oligomer_Left_1
        kf_TL1.name    = 'kf_TL1';
        kf_TL2.value   = kf_3d;     %3D Template_Right + Oligomer_Left_2
        kf_TL2.name    = 'kf_TL2';

        kb_TR1.value   = 71600;  %2D Template_Left + Oligomer_Right_1
        kb_TR1.name    = 'kb_TR1';
        kb_TR2.value   = 158;  %3D Template_Left + Oligomer_Right_2
        kb_TR2.name    = 'kb_TR2';
        kb_TL1.value   = 68000;    %2D Template_Right + Oligomer_Left_1
        kb_TL1.name    = 'kb_TL1';
        kb_TL2.value   = 75.5;    %3D Template_Right + Oligomer_Left_2
        kb_TL2.name    = 'kb_TL2';


        kf_lig.value     = k_lig;    % Nicked_Right -> Duplex
        kf_lig.name      = 'kf_lig';
        kf_duplex.value  = kf_2d;    % Simplexes -> Duplex
        kf_duplex.name   = 'kf_duplex';
        kb_duplex.value  = 0.37;  % Duplex -> Simplexes
        kb_duplex.name   = 'kb_duplex';

        % Background reaction constants
        kf_O1.value  = kf_2d;      %2D R1 + L1 -> O1
        kf_O1.name   = 'kf_O1';
        kf_O2.value  = 2e7;       %3D R2 + L2 -> O2 
        kf_O2.name   = 'kf_O2';
        kb_O1.value  = 56000;     %2D O1 -> R1 + L1
        kb_O1.name   = 'kb_O1';
        kb_O2.value  = 163.8;      %3D O2 -> R2 + L2 
        kb_O2.name   = 'kb_O2';

        % not happening in a twodimensional LIDA
        % kpara = 2e-2;                                 % Ligation of the two O1+- and O2+- 

        %dummy rates for reversible
        kb_lig.value  = 0;
        kb_lig.name   = 'kb_lig';

parameters = [kf_TR1,kf_TR2,kf_TL1,kf_TL2,kb_TR1,kb_TR2,kb_TL1,kb_TL2,kf_lig,...
              kf_duplex,kb_duplex,kf_O1,kf_O2,kb_O1,kb_O2,kb_lig];

end


%% ===================================================  fine exploration
if kset == "explore"

k_lig = 2e-2;
kf_2d = 1.5e8;

        kf_TR1.value   = kf_2d;    %2D Template_Left + Oligomer_Right_1
        kf_TR1.name    = 'kf_TR1';
        kf_TR2.value   = 2e7;     %3D Template_Left + Oligomer_Right_2
        kf_TR2.name    = 'kf_TR2';
        kf_TL1.value   = kf_2d;    %2D Template_Right + Oligomer_Left_1
        kf_TL1.name    = 'kf_TL1';
        kf_TL2.value   = 2e7;     %3D Template_Right + Oligomer_Left_2
        kf_TL2.name    = 'kf_TL2';

        kb_TR1.value   = 71600;  %2D Template_Left + Oligomer_Right_1
        kb_TR1.name    = 'kb_TR1';
        kb_TR2.value   = 158;  %3D Template_Left + Oligomer_Right_2
        kb_TR2.name    = 'kb_TR2';
        kb_TL1.value   = 68000;    %2D Template_Right + Oligomer_Left_1
        kb_TL1.name    = 'kb_TL1';
        kb_TL2.value   = 75.5;    %3D Template_Right + Oligomer_Left_2
        kb_TL2.name    = 'kb_TL2';


        kf_lig.value     = k_lig;    % Nicked_Right -> Duplex
        kf_lig.name      = 'kf_lig';
        kf_duplex.value  = kf_2d;    % Simplexes -> Duplex
        kf_duplex.name   = 'kf_duplex';
        kb_duplex.value  = 0.37;  % Duplex -> Simplexes
        kb_duplex.name   = 'kb_duplex';

        % Background reaction constants
        kf_O1.value  = kf_2d;      %2D R1 + L1 -> O1
        kf_O1.name   = 'kf_O1';
        kf_O2.value  = 2e7;       %3D R2 + L2 -> O2 
        kf_O2.name   = 'kf_O2';
        kb_O1.value  = 56000;     %2D O1 -> R1 + L1
        kb_O1.name   = 'kb_O1';
        kb_O2.value  = 163.8;      %3D O2 -> R2 + L2 
        kb_O2.name   = 'kb_O2';

        % not happening in a twodimensional LIDA
        % kpara = 2e-2;                                 % Ligation of the two O1+- and O2+- 

        %dummy rates for reversible
        kb_lig.value  = 0;
        kb_lig.name   = 'kb_lig';

parameters = [kf_TR1,kf_TR2,kf_TL1,kf_TL2,kb_TR1,kb_TR2,kb_TL1,kb_TL2,kf_lig,...
              kf_duplex,kb_duplex,kf_O1,kf_O2,kb_O1,kb_O2,kb_lig];

end

%% ===================================================  coarse grained exploration - freekbduplex

if kset == "mix-constrained"

kf_2d.value         = 2e7;
kf_2d.name          = 'kf_2d';
kf_mix.value        = 2e7;
kf_mix.name         = 'kf_mix';
kf_3d.value         = 2e7;
kf_3d.name          = 'kf_3d';

%defect asimmetry factor
defect_asimmetry    = 5e3;

kb_oligo_surf.value = 6e5;
kb_oligo_surf.name  = 'kb_oligo_surf';

kb_oligo_bulk.value = kb_oligo_surf.value / defect_asimmetry;
kb_oligo_bulk.name  = 'kb_oligo_bulk';

correction_initiation = 25;          %e^(-1.9/RT)
correction_skipped    = 20;          %anywhere from 5 to 67 which correspond to 
                                     %-1.0 to -2.5 kcal/mol free energy
ke_oligo_surf = kf_2d.value / kb_oligo_surf.value
ke_oligo_bulk = kf_mix.value / kb_oligo_bulk.value
ke_duplex = ke_oligo_surf * ke_oligo_bulk * correction_initiation * correction_skipped;
kb_duplex.value = kf_2d.value / ke_duplex
kb_duplex.name  = 'kb_duplex';

kf_lig.value  = 2e-2;
kf_lig.name   = 'kf_lig';
kb_lig.value  = 0;
kb_lig.name   = 'kb_lig';

parameters = [kf_2d,kf_3d,kf_mix,kb_duplex,kb_oligo_surf,kb_oligo_bulk,kf_lig];
 
end


%% ===================================================  coarse grained exploration - freekbduplex

if kset == "mix-free"

kf_2d.value         = 2e7;
kf_2d.name          = 'kf_2d';
kf_2d_duplex.value  = 2e7;
kf_2d_duplex.name   = 'kf_2d_duplex';
kf_mix.value        = 2e7;
kf_mix.name         = 'kf_mix';
kf_3d.value         = 2e7;
kf_3d.name          = 'kf_3d';

% ke_duplex       = 2e7;
% kb_duplex.value = kf_2d.value/ke_duplex;
kb_duplex.value = 0.3;
kb_duplex.name  = 'kb_duplex';

% kz_od             = 100;
kb_oligo_surf.value = 6e5;
kb_oligo_surf.name  = 'kb_oligo_surf';

defect_asimmetry    = 5e3;
kb_oligo_bulk.value = defect_asimmetry * kb_oligo_surf.value;
kb_oligo_bulk.name  = 'kb_oligo_bulk';

kf_lig.value        = 2e-2;
kf_lig.name         = 'kf_lig';


parameters = [kf_2d,kf_2d_duplex,kf_3d,kf_mix,kb_duplex,kb_oligo_surf,kb_oligo_bulk,kf_lig];

%         kf_TR1.value   = kf_2d.value;    %2D Template_Left + Oligomer_Right_1
%         kf_TR1.name    = 'kf_TR1';
%         kf_TR2.value   = kf_mix.value;     %3D Template_Left + Oligomer_Right_2
%         kf_TR2.name    = 'kf_TR2';
%         kf_TL1.value   = kf_2d.value;    %2D Template_Right + Oligomer_Left_1
%         kf_TL1.name    = 'kf_TL1';
%         kf_TL2.value   = kf_mix.value;     %3D Template_Right + Oligomer_Left_2
%         kf_TL2.name    = 'kf_TL2';
% 
%         kb_TR1.value   = 2e5;  %2D Template_Left + Oligomer_Right_1
%         kb_TR1.name    = 'kb_TR1';
%         kb_TR2.value   = 2e5;  %3D Template_Left + Oligomer_Right_2
%         kb_TR2.name    = 'kb_TR2';
%         kb_TL1.value   = 2e5;    %2D Template_Right + Oligomer_Left_1
%         kb_TL1.name    = 'kb_TL1';
%         kb_TL2.value   = 2e5;    %3D Template_Right + Oligomer_Left_2
%         kb_TL2.name    = 'kb_TL2';
% 
% 
%         kf_lig.value       = k_lig;    % Nicked_Right -> Duplex
%         kf_lig.name        = 'kf_lig';
%         kf_duplex.value    = kf_2d.value;    % Simplexes -> Duplex
%         kf_duplex.name     = 'kf_duplex';
%         kb_duplex.value    = kf_duplex.value / kz_duplex.value;  % Duplex -> Simplexes
%         kb_duplex.name     = 'kb_duplex';
% 
%         % Background reaction constants
%         kf_O1.value  = kf_2d.value;      %2D R1 + L1 -> O1
%         kf_O1.name   = 'kf_O1';
%         kf_O2.value  = kf_3d.value;       %3D R2 + L2 -> O2 
%         kf_O2.name   = 'kf_O2';
%         kb_O1.value  = 5e3;     %2D O1 -> R1 + L1
%         kb_O1.name   = 'kb_O1';
%         kb_O2.value  = 5e3;      %3D O2 -> R2 + L2 
%         kb_O2.name   = 'kb_O2';
% 
%         % not happening in a twodimensional LIDA
%         % kpara = 2e-2;                                 % Ligation of the two O1+- and O2+- 
% 
        %dummy rates for reversible
        kb_lig.value  = 0;
        kb_lig.name   = 'kb_lig';


end



