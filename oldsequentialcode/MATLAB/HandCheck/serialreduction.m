clear all


model = sbiomodel('parallel');
addcompartment(model,'bulk');

kf1.value = 1e3;
kf1.name = 'kf1';
kb1.value = 1e3;
kb1.name = 'kb1';

kf2.value = 2e5;
kf2.name = 'kf2';
kb2.value = 1e4;
kb2.name = 'kb2';

addparameter(model, kf1.name, kf1.value, 'units','1/second');
addparameter(model, kb1.name, kb1.value, 'units','1/second');
addparameter(model, kf2.name, kf2.value, 'units','1/second');
addparameter(model, kb2.name, kb2.value, 'units','1/second');

reaction1 = addreaction(model, "A <-> B");
kinlaw1 = addkineticlaw(reaction1, "MassAction");
set(kinlaw1, 'parametervariablenames',{kf1.name kb1.name});

reaction2 = addreaction(model, "B <-> C");
kinlaw2 = addkineticlaw(reaction2, "MassAction");
set(kinlaw2, 'parametervariablenames',{kf2.name kb2.name});

redmodel = sbiomodel('reducted');
addcompartment(redmodel,'bulk');

% kf_red.value = kf1.value - kb1.value*kf1.value/(kb1.value + kf2.value);
% kf_red.name = 'kf_red';
% kb_red.value = kb1.value*kb2.value/(kb1.value+kf2.value) + 20;
% kb_red.name = 'kb_red';

kf_red.value = kf1.value*kf2.value/(kb1.value + kf2.value);
kf_red.name = 'kf_red';
kb_red.value = kb1.value*kb2.value/(kb1.value + kf2.value);
kb_red.name = 'kb_red';

addparameter(redmodel, kf_red.name, kf_red.value, 'units','1/second');
addparameter(redmodel, kb_red.name, kb_red.value, 'units','1/second');

reducedreaction = addreaction(redmodel, "A <-> C");
reducedkinlaw = addkineticlaw(reducedreaction, "MassAction");
set(reducedkinlaw,'parametervariablenames',{kf_red.name kb_red.name});

model.species(1).initialamount = 1;
redmodel.species(1).initialamount = 1;

abstol = 1e-20;
reltol = 1e-12;
timeunits = 'second';
runtime = 1e-2;

conf = model.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(model);

conf = redmodel.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(redmodel);

sbioplot(sbiosimulate(model));
hold on
sbioplot(sbiosimulate(redmodel));

