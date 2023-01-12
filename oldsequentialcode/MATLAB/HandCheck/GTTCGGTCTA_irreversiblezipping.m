% GTTCGGTCTA - main strand
% CAAGCCAGAT - complementary
clear all

%% Initial Model Specification 

kf_exp = 1.14594998e6;

%my python script spits out this value for the geometric rate
kf_dl  = 1.676959e7*0.7;
nonspecific_correction = 1;

% from porschke1973
kf_zipping = 2e9;
sliding_correction = 1;

dangling_correction = 1e-1;

model = sbiomodel('GTT');
bulk = addcompartment(model,'bulk');

nodes_names = ["N1", "N2", "N3", "N4", "N5", "N6", "N7"];

%% Reaction 1
reaction_a = addreaction(model, "SS + SS <-> N1");
kin_a = addkineticlaw(reaction_a, 'MassAction');
kf_a.value = kf_dl*nonspecific_correction;
kf_a.name = 'kf_a';
kb_a.value = 9.351890e+05*dangling_correction; 
kb_a.name = 'kb_a';
addparameter(model, kf_a.name, kf_a.value, 'units', '1/(molarity*second)');
addparameter(model, kb_a.name, kb_a.value, 'units', '1/second');
set(kin_a, 'parametervariablenames',{kf_a.name kb_a.name});

%% Reaction 2
reaction_b = addreaction(model, "SS + SS <-> N2");
kin_b = addkineticlaw(reaction_b, 'MassAction');
kf_b.value = kf_dl*nonspecific_correction;
kf_b.name = 'kf_b';
kb_b.value = 3.250175e+05;
kb_b.name = 'kb_b';
addparameter(model, kf_b.name, kf_b.value, 'units', '1/(molarity*second)');
addparameter(model, kb_b.name, kb_b.value, 'units', '1/second');
set(kin_b, 'parametervariablenames',{kf_b.name kb_b.name});

%% Reaction 3
reaction_c = addreaction(model, "SS + SS <-> N3");
kin_c = addkineticlaw(reaction_c, 'MassAction');
kf_c.value = kf_dl*nonspecific_correction;
kf_c.name = 'kf_c';
kb_c.value = 8.197089e+04*dangling_correction;
kb_c.name = 'kb_c';
addparameter(model, kf_c.name, kf_c.value, 'units', '1/(molarity*second)');
addparameter(model, kb_c.name, kb_c.value, 'units', '1/second');
set(kin_c, 'parametervariablenames',{kf_c.name kb_c.name});

%% Reaction 4
reaction_d = addreaction(model, "SS + SS <-> N4");
kin_d = addkineticlaw(reaction_d, 'MassAction');
kf_d.value = kf_dl*nonspecific_correction;
kf_d.name = 'kf_d';
kb_d.value = 6.043337e+04*dangling_correction; 
kb_d.name = 'kb_d';
addparameter(model, kf_d.name, kf_d.value, 'units', '1/(molarity*second)');
addparameter(model, kb_d.name, kb_d.value, 'units', '1/second');
set(kin_d, 'parametervariablenames',{kf_d.name kb_d.name});

%% Reaction 5
reaction_e = addreaction(model, "SS + SS <-> N5");
kin_e = addkineticlaw(reaction_e, 'MassAction');
kf_e.value = kf_dl*nonspecific_correction;
kf_e.name = 'kf_e';
kb_e.value = 2.358589e+05*dangling_correction;
kb_e.name = 'kb_e';
addparameter(model, kf_e.name, kf_e.value, 'units', '1/(molarity*second)');
addparameter(model, kb_e.name, kb_e.value, 'units', '1/second');
set(kin_e, 'parametervariablenames',{kf_e.name kb_e.name});

%% Reaction 6
reaction_f = addreaction(model, "SS + SS <-> N6");
kin_f = addkineticlaw(reaction_f, 'MassAction');
kf_f.value = kf_dl*nonspecific_correction;
kf_f.name = 'kf_f';
kb_f.value = 7.218884e+05*dangling_correction;
kb_f.name = 'kb_f';
addparameter(model, kf_f.name, kf_f.value, 'units', '1/(molarity*second)');
addparameter(model, kb_f.name, kb_f.value, 'units', '1/second');
set(kin_f, 'parametervariablenames',{kf_f.name kb_f.name});

%% Reaction 7 
reaction_g = addreaction(model, "SS + SS <-> N7");
kin_g = addkineticlaw(reaction_g, 'MassAction');
kf_g.value = kf_dl*nonspecific_correction; 
kf_g.name = 'kf_g';
kb_g.value = 3.927066e+06*dangling_correction;
kb_g.name = 'kb_g';
addparameter(model, kf_g.name, kf_g.value, 'units', '1/(molarity*second)');
addparameter(model, kb_g.name, kb_g.value, 'units', '1/second'); 
set(kin_g, 'parametervariablenames',{kf_g.name kb_g.name});


%% Zipping General
kf_z1.value = kf_zipping; 
kf_z1.name = 'kf_z1';
addparameter(model, kf_z1.name, kf_z1.value, 'units', '1/second');

%% Zipping 1
reaction_z1 = addreaction(model, "N1 -> D");
kin_z1 = addkineticlaw(reaction_z1, 'MassAction');
set(kin_z1, 'parametervariablenames',{kf_z1.name});
%% Zipping 2
reaction_z2 = addreaction(model, "N2 -> D");
kin_z2 = addkineticlaw(reaction_z2, 'MassAction');
set(kin_z2, 'parametervariablenames',{kf_z1.name});
%% Zipping 3
reaction_z3 = addreaction(model, "N3 -> D");
kin_z3 = addkineticlaw(reaction_z3, 'MassAction');
set(kin_z3, 'parametervariablenames',{kf_z1.name});
%% Zipping 4
reaction_z4 = addreaction(model, "N4 -> D");
kin_z4 = addkineticlaw(reaction_z4, 'MassAction');
set(kin_z4, 'parametervariablenames',{kf_z1.name});
%% Zipping 5
reaction_z5 = addreaction(model, "N5 -> D");
kin_z5 = addkineticlaw(reaction_z5, 'MassAction');
set(kin_z5, 'parametervariablenames',{kf_z1.name});
%% Zipping 6
reaction_z6 = addreaction(model, "N6 -> D");
kin_z6 = addkineticlaw(reaction_z6, 'MassAction');
set(kin_z6, 'parametervariablenames',{kf_z1.name});
%% Zipping 7
reaction_z7 = addreaction(model, "N7 -> D");
kin_z7 = addkineticlaw(reaction_z7, 'MassAction');
set(kin_z7, 'parametervariablenames',{kf_z1.name});

%% General Backward Reaction
reaction_back = addreaction(model, "D -> SS + SS");
kin_back = addkineticlaw(reaction_back, "MassAction");
kf_back.name = 'kf_back';
kf_back.value = 1;
addparameter(model, kf_back.name, kf_back.value, "Units", "1/second");
set(kin_back, 'parametervariablenames',{kf_back.name});

%% Model Execution 

% Initial SS concentration 
ss0 = 1e-2;

model.species(1).initialamount = ss0;

abstol = 1e-20;
reltol = 1e-12;
timeunits = 'second';
runtime = 5e-5;

conf = model.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(model);

sim_GTT = sbiosimulate(model);
time = sim_GTT.time;
duplex = sim_GTT.Data(:,9);


fitmodel = sbiomodel('twostate');
reaction = addreaction(fitmodel,'SS + SS <-> D');
set(fitmodel.compartment,'units','liter')
initial = ss0;
fitmodel.species(1).initialamount = initial;
set(fitmodel.species, 'units', 'molarity');
kinlaw = addkineticlaw(reaction,'MassAction');
addparameter(fitmodel, 'kf', 1, 'units', '1/(molarity*second)');
addparameter(fitmodel, 'kb', 1, 'units', '1/second');
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