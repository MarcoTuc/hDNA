% ACACACACCA - main strand
% TGTGTGTGGT - complementary
clear all


kf_exp = 4.15934656e6;

%my python script spits out this value for the geometric rate
kf_dl  = 1.676959e7;

% from porschke1973
kf_sliding = 5e3;

model = sbiomodel('ACA');
bulk = addcompartment(model,'bulk');

nodes_names = ["SS", "L1", "L2", "R1", "D"];


reaction_a = addreaction(model, "SS + SS <-> L1");
kin_a = addkineticlaw(reaction_a, 'MassAction');
kf_a.value = kf_dl;
kf_a.name = 'kf_a';
kb_a.value = 3.427038e+05; 
kb_a.name = 'kb_a';
addparameter(model, kf_a.name, kf_a.value, 'units', '1/(molarity*second)');
addparameter(model, kb_a.name, kb_a.value, 'units', '1/second');
set(kin_a, 'parametervariablenames',{kf_a.name kb_a.name});


reaction_b = addreaction(model, "L1 <-> L2");
kin_b = addkineticlaw(reaction_b, 'MassAction');
kf_b.value = kf_sliding;
kf_b.name = 'kf_b';
kb_b.value = 1.589447e+01;
kb_b.name = 'kb_b';
addparameter(model, kf_b.name, kf_b.value, 'units', '1/second');
addparameter(model, kb_b.name, kb_b.value, 'units', '1/second');
set(kin_b, 'parametervariablenames',{kf_b.name kb_b.name});


reaction_c = addreaction(model, "L2 <-> D");
kin_c = addkineticlaw(reaction_c, 'MassAction');
kf_c.value = 5e3;
kf_c.name = 'kf_c';
kb_c.value = 9.780507e-03;
kb_c.name = 'kb_c';
addparameter(model, kf_c.name, kf_c.value, 'units', '1/second');
addparameter(model, kb_c.name, kb_c.value, 'units', '1/second');
set(kin_c, 'parametervariablenames',{kf_c.name kb_c.name});


reaction_d = addreaction(model, "SS + SS <-> D");
kin_d = addkineticlaw(reaction_d, 'MassAction');
kf_d.value = kf_dl;
kf_d.name = 'kf_d';
kb_d.value = 2.131013e-03; 
kb_d.name = 'kb_d';
addparameter(model, kf_d.name, kf_d.value, 'units', '1/(molarity*second)');
addparameter(model, kb_d.name, kb_d.value, 'units', '1/second');
set(kin_d, 'parametervariablenames',{kf_d.name kb_d.name});


reaction_e = addreaction(model, "SS + SS <-> L2");
kin_e = addkineticlaw(reaction_e, 'MassAction');
kf_e.value = kf_dl;
kf_e.name = 'kf_e';
kb_e.value = 1.089419e+03;
kb_e.name = 'kb_e';
addparameter(model, kf_e.name, kf_e.value, 'units', '1/(molarity*second)');
addparameter(model, kb_e.name, kb_e.value, 'units', '1/second');
set(kin_e, 'parametervariablenames',{kf_e.name kb_e.name});


reaction_f = addreaction(model, "SS + SS <-> R1");
kin_f = addkineticlaw(reaction_f, 'MassAction');
kf_f.value = kf_dl;
kf_f.name = 'kf_f';
kb_f.value = 5.237242e+03;
kb_f.name = 'kb_f';
addparameter(model, kf_f.name, kf_f.value, 'units', '1/(molarity*second)');
addparameter(model, kb_f.name, kb_f.value, 'units', '1/second');
set(kin_f, 'parametervariablenames',{kf_f.name kb_f.name});


reaction_g = addreaction(model, "R1 <-> D");
kin_g = addkineticlaw(reaction_g, 'MassAction');
kf_g.value = kf_sliding; 
kf_g.name = 'kf_g';
kb_g.value = 2.034481e-03;
kb_g.name = 'kb_g';
addparameter(model, kf_g.name, kf_g.value, 'units', '1/second');
addparameter(model, kb_g.name, kb_g.value, 'units', '1/second');
set(kin_g, 'parametervariablenames',{kf_g.name kb_g.name});

ss0 = 1;

model.species(1).initialamount = ss0;

abstol = 1e-20;
reltol = 1e-12;
timeunits = 'second';
runtime = 1e-3;

conf = model.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(model);

sim_ACA = sbiosimulate(model);
time = sim_ACA.time;
duplex = sim_ACA.Data(:,4);


fitmodel = sbiomodel('twostate');
reaction = addreaction(fitmodel,'SS + SS <-> D');
set(fitmodel.compartment,'units','liter')
initial = 1;
fitmodel.species(1).initialamount = initial;
set(fitmodel.species, 'units', 'molarity');
kinlaw = addkineticlaw(reaction,'MassAction');
addparameter(fitmodel, 'kf', 1, 'units', '1/(molarity*second)');
addparameter(fitmodel, 'kb', 1, 'units', '1/second');
set(kinlaw,'parametervariablenames',{'kf' 'kb'});

G = groupedData(table(time,duplex));
responsemap = 'D = duplex';
estimate = estimatedInfo('kf');
fitresult = sbiofit(fitmodel, G, responsemap, estimate);

[yfit,paramEstim] = fitted(fitresult);


sbioplot(sim_ACA);
hold  on
sbioplot(yfit);

fitresult.parameterestimate