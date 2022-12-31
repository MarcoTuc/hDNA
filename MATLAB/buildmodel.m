
%#ok<*AGROW>
% params_formulas

function model = buildmodel(name,kset)

model = sbiomodel(name);
membrane = addcompartment(model, 'membrane');

rate_constants;


%% Species routine
if kset == "explore-stoch"
species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];
y0_olig         = 100;
y0_temp         = 100;
species_initial = [y0_temp, y0_olig, y0_olig, 0, 0, 0, 0, y0_olig, y0_olig, 0, 0, 0, 0, 0, 0];
species_initial = double(species_initial);
else
species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];
y0_olig         = 1e-4;
y0_temp         = 1e-9;
species_initial = [y0_temp, y0_olig, y0_olig, 0, 0, 0, 0, y0_olig, y0_olig, 0, 0, 0, 0, 0, 0];
species_initial = double(species_initial);
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
    addreaction(model,reactionstring);
end


%% Rates routine

if kset == "explore" || kset == "realistic" || kset == "explore-stoch"
    fwd_rates = [kf_TR1,kf_TR2,kf_TR2,kf_TR1,kf_lig,kb_duplex,kf_TL1,kf_TL2,kf_TL2,kf_TL1,kf_lig,kf_O1,kf_O2];
    bwd_rates = [kb_TR1,kb_TR2,kb_TR2,kb_TR1,kb_lig,kf_duplex,kb_TL1,kb_TL2,kb_TL2,kb_TL1,kb_lig,kb_O1,kb_O2];
end

if kset == "mix-free"
    fwd_rates = [kf_2d,kf_mix,kf_mix,kf_2d,kf_lig,kb_duplex,kf_2d,kf_mix,kf_mix,kf_2d,kf_lig,kf_2d,kf_3d];
    bwd_rates = [kb_oligo_surf,kb_oligo_bulk,kb_oligo_bulk,kb_oligo_surf,kb_lig,kf_2d_duplex,kb_oligo_surf,kb_oligo_bulk,kb_oligo_bulk,kb_oligo_surf,kb_lig,kb_oligo_surf,kb_oligo_bulk];
end

if kset == "mix-constrained"
    fwd_rates = [kf_2d,         kf_mix,        kf_mix,        kf_2d,         kf_lig, kb_duplex, kf_2d,         kf_mix,        kf_mix,        kf_2d,         kf_lig, kf_2d,         kf_3d];
    bwd_rates = [kb_oligo_surf, kb_oligo_bulk, kb_oligo_bulk, kb_oligo_surf, kb_lig, kf_2d,     kb_oligo_surf, kb_oligo_bulk, kb_oligo_bulk, kb_oligo_surf, kb_lig, kb_oligo_surf, kb_oligo_bulk];
end

fwd_rates_matrix = cell2mat({fwd_rates.value});
bwd_rates_matrix = cell2mat({bwd_rates.value});

F = Stochio*diag(fwd_rates_matrix);
B = Stochio*diag(fwd_rates_matrix);

kineticmodel = 'MassAction';

for p = 1:length(parameters)
    addparameter(model, parameters(p).name, parameters(p).value);
end

for r = 1:length(model.reactions)
    kin = addkineticlaw(model.reactions(r),kineticmodel);
    if model.reaction(r).reversible
        set(kin,'parametervariablenames',{fwd_rates(r).name bwd_rates(r).name});
    else
        set(model.reactions(r).kineticlaw,'parametervariablenames',{fwd_rates(r).name});
    end
end

v1 = sbiovariant('v1');
addvariant(model,'v1');

cursorvar = sbiovariant('cursorvar');
addvariant(model,'cursorvar');

clearvars -except model v1 cursorvar Stochio species_names

end