% GTTCGGTCTA - main strand
% CAAGCCAGAT - complementary

%% Initial Model Specification 

clear all

kf_exp = 1.14e6;

%my python script spits out this value for the geometric rate

model = sbiomodel('GTT');
bulk = addcompartment(model,'bulk');

% Reaction Nucleation
reaction_a = addreaction(model, "SS + SS <-> N1");
kin_a = addkineticlaw(reaction_a, 'MassAction');
nonspecific_correction = 1;
kf_a.value = 2e7*nonspecific_correction;
kf_a.name = 'kf_a';
kb_a.value = 2e7; 
kb_a.name = 'kb_a';
addparameter(model, kf_a.name, kf_a.value, 'units', '1/(molarity*second)');
addparameter(model, kb_a.name, kb_a.value, 'units', '1/second');
set(kin_a, 'parametervariablenames',{kf_a.name kb_a.name});

% Reaction Zipping
reaction_z = addreaction(model, "N1 <-> D");
kin_z = addkineticlaw(reaction_z, 'MassAction');
kf_z.value = 2e9;
kf_z.name = 'kf_z';
kb_z.value = 4e8;
kb_z.name = 'kb_z';
addparameter(model, kf_z.name, kf_z.value, 'units', '1/second');
addparameter(model, kb_z.name, kb_z.value, 'units', '1/second');
set(kin_z, 'parametervariablenames',{kf_z.name kb_z.name});


%% Model Execution 

% Initial SS concentration 
ss0 = 1e-2;

model.species(1).initialamount = ss0;

abstol = 1e-20;
reltol = 1e-12;
timeunits = 'second';
runtime = 1e-5;

conf = model.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(model);

sim_GTT = sbiosimulate(model);
time = sim_GTT.time;
duplex = sim_GTT.Data(:,3);

kb = 1;

fitmodel = sbiomodel('twostate');
reaction = addreaction(fitmodel,'SS + SS <-> D');
set(fitmodel.compartment,'units','liter')
initial = ss0;
fitmodel.species(1).initialamount = initial;
set(fitmodel.species, 'units', 'molarity');
kinlaw = addkineticlaw(reaction,'MassAction');
addparameter(fitmodel, 'kf', 1, 'units', '1/(molarity*second)');
addparameter(fitmodel, 'kb', kb, 'units', '1/second');
set(kinlaw,'parametervariablenames',{'kf' 'kb'});

G = groupedData(table(time,duplex));
responsemap = 'D = duplex';
estimate = estimatedInfo('kf');
fitresult = sbiofit(fitmodel, G, responsemap, estimate, [], 'nlinfit');

[yfit,paramEstim] = fitted(fitresult);


sbioplot(sim_GTT);
hold  on
sbioplot(yfit);

fitresult.ParameterEstimates