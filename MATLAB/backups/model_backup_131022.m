

%#ok<*AGROW>
params_formulas

lida = sbiomodel('LIDA');
membrane = addcompartment(lida, 'membrane');

k_set = "explore";
rate_constants;

stochastic = "false";

max_data = []

%% Species routine


species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];
y0_olig         = 1e-4;
y0_temp         = 1e-9;
species_initial = [y0_temp, y0_olig, y0_olig, 0, 0, 0, 0, y0_olig, y0_olig, 0, 0, 0, 0, 0, 0];
% species_initial = [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0];
% species_initial = rand(1,15)
species_initial = double(species_initial);

    if stochastic == "true"
        species_initial = species_initial*NA
    end


for i = 1:length(species_names)
species = addspecies(membrane,species_names(i),species_initial(i));
end

%% Reactions routine
Rev_matrix = [
      1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1];
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


for r = 1:size(Stochio,2)
    reactants = [];
    products = [];
        for i = 1:max(size(Stochio))
            if Stochio(i,r) == -1
               reactants = [reactants species_names(i)];
            end
            if Stochio(i,r) == 1
               products = [products species_names(i)];
            end
        end
    if Rev_matrix(r) == 0
        reactionstring = strjoin(reactants,' + ')+' -> '+strjoin(products,' + ');
    else
        reactionstring = strjoin(reactants,' + ')+' <-> '+strjoin(products,' + ');
    end 
    addreaction(lida,reactionstring);
end


%% Rates routine




fwd_rates = [kf_TR1,kf_TR2,kf_TR2,kf_TR1,kf_lig_left,kb_duplex,kf_TL1,kf_TL2,kf_TL2,kf_TL1,kf_lig_right,kf_O1,kf_O2];
bwd_rates = [kb_TR1,kb_TR2,kb_TR2,kb_TR1,kb_lig_left,kf_duplex,kb_TL1,kb_TL2,kb_TL2,kb_TL1,kb_lig_right,kb_O1,kb_O2];

fwd_rates_matrix = cell2mat({fwd_rates.value});
bwd_rates_matrix = cell2mat({bwd_rates.value});

if stochastic == "true"
    fwd_rates_matrix = fwd_rates_matrix/NA
    bwd_rates_matrix = bwd_rates_matrix/NA
end

F = Stochio*diag(fwd_rates_matrix);
B = Stochio*diag(fwd_rates_matrix);

kineticmodel = 'MassAction';

for r = 1:length(lida.reactions)
    kin = addkineticlaw(lida.reactions(r),kineticmodel);
    if lida.reaction(r).reversible
        addparameter(kin,fwd_rates(r).name, fwd_rates(r).value);
        addparameter(kin,bwd_rates(r).name, bwd_rates(r).value);
        set(kin,'parametervariablenames',{fwd_rates(r).name bwd_rates(r).name});
    else
        addparameter(lida.reactions(r).kineticlaw,fwd_rates(r).name,  fwd_rates(r).value);
        set(lida.reactions(r).kineticlaw,'parametervariablenames',{fwd_rates(r).name});
    end
end



%% Simulate

abstol = 1e-20;
reltol = 1e-13;
runtime = 5e6;
timeunits = 'seconds';

conf = lida.config;

if stochastic ~= "true"
    set(conf, 'solvertype','ode15s');
    set(conf.solveroptions, 'absolutetolerance',abstol);
    set(conf.solveroptions, 'relativetolerance',reltol);
end
if stochastic == "true"
    set(conf, 'solvertype','ode15s');

end

% set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(lida);
SD = sbiosimulate(lida);

sbioplot(SD)

%% Find Max of Nl and Nr

% Nl = SD.Data(:,6);
% Nr = SD.Data(:,12);
% [max_Nl,indx_Nl] = max(Nl);
% [max_Nr,indx_Nr] = max(Nr);
% tmax_Nl = SD.Time(indx_Nl);
% tmax_Nr = SD.Time(indx_Nr);
% 
% figure(5)
% plot(SD.Time,Nl)
% hold on
% plot(SD.Time,Nr)
% hold on 
% plot(tmax_Nl,max_Nl,'r+')
% plot(tmax_Nr,max_Nr,'b+')

%% max data 

[max_N, indx_N] = max(SD.Data(:,6));
tmax_N  = SD.Time(indx_N);
new_data = [max_N tmax_N];
max_data = [max_data new_data];

kf_d = lida.reactions(6).kineticlaw.parameters(2).value;
kb_d = lida.reactions(6).kineticlaw.parameters(1).value;
k_lig = lida.reactions(5).kineticlaw.parameters(1).value;
nl = SD.Data(end,5);
nr = SD.Data(end,12);
Tl = SD.Data(end,1);
Tr = SD.Data(end,7);
D = SD.Data(end,13);

KED = (kf_d/kb_d) + (k_lig/kb_d)*((nl+nr)/(Tl*Tr));
KEDDDD = D / (Tl*Tr);
Kminus = kf_d + k_lig*((nl+nr)/(Tl*Tr));

tau = (k_lig*y0_olig);
% xline(tau)