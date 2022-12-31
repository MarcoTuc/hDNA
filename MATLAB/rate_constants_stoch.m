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

        kb_TR1.value   = 71600 / Na;  %2D Template_Left + Oligomer_Right_1
        kb_TR1.name    = 'kb_TR1';
        kb_TR2.value   = 158 / Na;  %3D Template_Left + Oligomer_Right_2
        kb_TR2.name    = 'kb_TR2';
        kb_TL1.value   = 68000 / Na;    %2D Template_Right + Oligomer_Left_1
        kb_TL1.name    = 'kb_TL1';
        kb_TL2.value   = 75.5 / Na;    %3D Template_Right + Oligomer_Left_2
        kb_TL2.name    = 'kb_TL2';


        kf_lig.value     = k_lig;    % Nicked_Right -> Duplex
        kf_lig.name      = 'kf_lig';
        kf_duplex.value  = kf_2d;    % Simplexes -> Duplex
        kf_duplex.name   = 'kf_duplex';
        kb_duplex.value  = 0.37 / Na;  % Duplex -> Simplexes
        kb_duplex.name   = 'kb_duplex';

        % Background reaction constants
        kf_O1.value  = kf_2d;      %2D R1 + L1 -> O1
        kf_O1.name   = 'kf_O1';
        kf_O2.value  = 2e7 / Na;       %3D R2 + L2 -> O2 
        kf_O2.name   = 'kf_O2';
        kb_O1.value  = 56000 / Na;     %2D O1 -> R1 + L1
        kb_O1.name   = 'kb_O1';
        kb_O2.value  = 163.8 / Na;      %3D O2 -> R2 + L2 
        kb_O2.name   = 'kb_O2';

        % not happening in a twodimensional LIDA
        % kpara = 2e-2;                                 % Ligation of the two O1+- and O2+- 

        %dummy rates for reversible
        kb_lig.value  = 0;
        kb_lig.name   = 'kb_lig';

parameters = [kf_TR1,kf_TR2,kf_TL1,kf_TL2,kb_TR1,kb_TR2,kb_TL1,kb_TL2,kf_lig,...
              kf_duplex,kb_duplex,kf_O1,kf_O2,kb_O1,kb_O2,kb_lig];

end