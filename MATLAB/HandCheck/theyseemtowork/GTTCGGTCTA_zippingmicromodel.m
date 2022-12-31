% GTTCGGTCTA - main strand
% CAAGCCAGAT - complementary

%% Initial Model Specification 

clear all

kf_exp = 1.14e6;

%my python script spits out this value for the geometric rate
kf_dl  = 4.214659e+06;
nonspecific_correction = 1/8;

dangling_correction = 1;

model = sbiomodel('GTT');
bulk = addcompartment(model,'bulk');


% Reaction 1
reaction_a = addreaction(model, "SS1 + SS2 <-> N1");
kin_a = addkineticlaw(reaction_a, 'MassAction');
kf_a.value = kf_dl*nonspecific_correction;
kf_a.name = 'kf_a';
kb_a.value = 1.315109e+05*dangling_correction; 
kb_a.name = 'kb_a';
addparameter(model, kf_a.name, kf_a.value, 'units', '1/(molarity*second)');
addparameter(model, kb_a.name, kb_a.value, 'units', '1/second');
set(kin_a, 'parametervariablenames',{kf_a.name kb_a.name});

% Reaction 2
reaction_b = addreaction(model, "SS1 + SS2 <-> N2");
kin_b = addkineticlaw(reaction_b, 'MassAction');
kf_b.value = kf_dl*nonspecific_correction;
kf_b.name = 'kf_b';
kb_b.value = 4.570558e+04*dangling_correction;
kb_b.name = 'kb_b';
addparameter(model, kf_b.name, kf_b.value, 'units', '1/(molarity*second)');
addparameter(model, kb_b.name, kb_b.value, 'units', '1/second');
set(kin_b, 'parametervariablenames',{kf_b.name kb_b.name});

% Reaction 3
reaction_c = addreaction(model, "SS1 + SS2 <-> N3");
kin_c = addkineticlaw(reaction_c, 'MassAction');
kf_c.value = kf_dl*nonspecific_correction;
kf_c.name = 'kf_c';
kb_c.value = 1.152716e+04*dangling_correction;
kb_c.name = 'kb_c';
addparameter(model, kf_c.name, kf_c.value, 'units', '1/(molarity*second)');
addparameter(model, kb_c.name, kb_c.value, 'units', '1/second');
set(kin_c, 'parametervariablenames',{kf_c.name kb_c.name});

% Reaction 4
reaction_d = addreaction(model, "SS1 + SS2 <-> N4");
kin_d = addkineticlaw(reaction_d, 'MassAction');
kf_d.value = kf_dl*nonspecific_correction;
kf_d.name = 'kf_d';
kb_d.value = 8.498443e+03*dangling_correction; 
kb_d.name = 'kb_d';
addparameter(model, kf_d.name, kf_d.value, 'units', '1/(molarity*second)');
addparameter(model, kb_d.name, kb_d.value, 'units', '1/second');
set(kin_d, 'parametervariablenames',{kf_d.name kb_d.name});

% Reaction 5
reaction_e = addreaction(model, "SS1 + SS2 <-> N5");
kin_e = addkineticlaw(reaction_e, 'MassAction');
kf_e.value = kf_dl*nonspecific_correction;
kf_e.name = 'kf_e';
kb_e.value = 3.316766e+04*dangling_correction;
kb_e.name = 'kb_e';
addparameter(model, kf_e.name, kf_e.value, 'units', '1/(molarity*second)');
addparameter(model, kb_e.name, kb_e.value, 'units', '1/second');
set(kin_e, 'parametervariablenames',{kf_e.name kb_e.name});

% Reaction 6
reaction_f = addreaction(model, "SS1 + SS2 <-> N6");
kin_f = addkineticlaw(reaction_f, 'MassAction');
kf_f.value = kf_dl*nonspecific_correction;
kf_f.name = 'kf_f';
kb_f.value = 1.015156e+05*dangling_correction;
kb_f.name = 'kb_f';
addparameter(model, kf_f.name, kf_f.value, 'units', '1/(molarity*second)');
addparameter(model, kb_f.name, kb_f.value, 'units', '1/second');
set(kin_f, 'parametervariablenames',{kf_f.name kb_f.name});

% Reaction 7 
reaction_g = addreaction(model, "SS1 + SS2 <-> N7");
kin_g = addkineticlaw(reaction_g, 'MassAction');
kf_g.value = kf_dl*nonspecific_correction; 
kf_g.name = 'kf_g';
kb_g.value = 5.522437e+05*dangling_correction;
kb_g.name = 'kb_g';
addparameter(model, kf_g.name, kf_g.value, 'units', '1/(molarity*second)');
addparameter(model, kb_g.name, kb_g.value, 'units', '1/second'); 
set(kin_g, 'parametervariablenames',{kf_g.name kb_g.name});


% Zipping Reactions

kf_newbp = 4e7;
activation_energy = 1.7; % Ea/kbT = 3 --> Ea = 3kbT (manghi destanville)

% kf_z.value = kf_newbp*exp(-activation_energy);
kf_z.value = kf_newbp;
kf_z.name = 'kf_z';
kb_z.value = 15e5;
kb_z.name = 'kb_z';

nbp = 10;
coresize = 3;

addparameter(model, kf_z.name,kf_z.value,"units","1/second");
addparameter(model, kb_z.name,kb_z.value,"units","1/second");

numcores = nbp - coresize;

for core = 1:numcores
cbase = coresize;
firstnucleus = coresize + 1;
reaction_string_fbp = string("N"+core+" <-> "+"N"+core+"+"+firstnucleus);
firstnewbp = addreaction(model, reaction_string_fbp);
kinlawfbp = addkineticlaw(firstnewbp,"MassAction");
set(kinlawfbp, 'ParameterVariableNames',{kf_z.name kb_z.name});
for base = coresize:nbp - 2
cbase = cbase+1;
if base == nbp - 2
reaction_string = string("N"+core+"+"+cbase+" <-> D");
else
nextbase = cbase+1;
reaction_string = string("N"+core+"+"+cbase+" <-> N"+core+"+"+nextbase);
end
reaction = addreaction(model, reaction_string);
kinlaw = addkineticlaw(reaction,"MassAction");
set(kinlaw,'ParameterVariableNames',{kf_z.name kb_z.name});
end
end


%% Model Execution 

% Initial SS concentration 
pseudoss0 = 1;
ss0 = 1;

model.species(1).initialamount = ss0;
model.species(2).initialamount = pseudoss0;

abstol = 1e-20;
reltol = 1e-12;
timeunits = 'second';
runtime = 1e-6;

conf = model.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(model);

sim_GTT = sbiosimulate(model);
time = sim_GTT.time;
duplex = sim_GTT.Data(:,16);

delay = 1;
time = time(delay:end);
delayedduplex = duplex(delay:end);

kb = 1e-4;

fitmodel = sbiomodel('twostate');
reaction = addreaction(fitmodel,'SS1 + SS2 <-> D');
set(fitmodel.compartment,'units','liter')
fitmodel.species(1).initialamount = ss0;
fitmodel.species(2).initialamount = pseudoss0;
set(fitmodel.species, 'units', 'molarity');
kinlaw = addkineticlaw(reaction,'MassAction');
addparameter(fitmodel, 'kf', 1, 'units', '1/(molarity*second)');
addparameter(fitmodel, 'kb', kb, 'units', '1/second');
set(kinlaw,'parametervariablenames',{'kf' 'kb'});

G = groupedData(table(time,delayedduplex));
responsemap = 'D = delayedduplex';
estimate = estimatedInfo('kf');
fitresult = sbiofit(fitmodel, G, responsemap, estimate, [], 'nlinfit');

[yfit,paramEstim] = fitted(fitresult);


sbioplot(sim_GTT);
hold  on
sbioplot(yfit);

fitresult.ParameterEstimates


%% See


plot(time,delayedduplex);