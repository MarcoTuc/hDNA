% GTTCGGTCTA - main strand
% CAAGCCAGAT - complementary

%% Initial Model Specification 
clear all

kf_newbp = 2e9;
activation_energy = 3; % Ea/kbT = 3 --> Ea = 3kbT (manghi destanville)

kf_z.value = kf_newbp*exp(-activation_energy);
kf_z.name = 'kf_z';
kb_z.value = 2e7;
kb_z.name = 'kb_z';

nbp = 10;
coresize = 3;

model = sbiomodel('ZIPPING');
bulk = addcompartment(model,'bulk');
addparameter(model, kf_z.name,kf_z.value,"units","1/second");
addparameter(model, kb_z.name,kb_z.value,"units","1/second");


zipnum = nbp - coresize;

for core = coresize:nbp-1
initialcore = core;
newcore = initialcore+1;
reaction_string = string("N"+initialcore+" <-> N"+newcore);
reaction = addreaction(model, reaction_string);
kinlaw = addkineticlaw(reaction,"MassAction");
set(kinlaw,'ParameterVariableNames',{kf_z.name kb_z.name});
end



%% Model Execution 

% Initial SS concentration 
ss0 = 1e-2;

model.species(1).initialamount = ss0;

abstol = 1e-20;
reltol = 1e-12;
timeunits = 'second';
runtime = 5e-7;

conf = model.config;
set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);
set(conf, 'stoptime',runtime);
sbioaccelerate(model);

sim_GTT = sbiosimulate(model);
time = sim_GTT.time;
duplex = sim_GTT.Data(:,end);

kb = 1;

fitmodel = sbiomodel('twostate');
reaction = addreaction(fitmodel,'SS + SS <-> D');
set(fitmodel.compartment,'units','liter')
initial = 2*ss0;
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