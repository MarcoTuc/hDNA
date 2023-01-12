kset = "explore-stoch";
modelname = "lida_stoch";
lida_stoch = buildmodel(modelname,kset)

species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];

%% Simulation Configuration 

abstol = 1e-14;
reltol = 1e-13;
runtime = 2e26;
timeunits = 'second';

conf = lida_stoch.config;

set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);

set(conf, 'stoptime',runtime);
sbioaccelerate(lida_stoch);
simdata = sbiosimulate(lida_stoch);
sbioplot(simdata);



