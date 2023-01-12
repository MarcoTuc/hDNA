%% Custom
clear all
m = sbiomodel('model');

addspecies(m,'A',1);
addspecies(m,'B',1);
addspecies(m,'C',0);

addreaction(m,'A + B <-> C');

kf = 100;
kz = 100;

delete(sbioselect('Type','kineticlaw','KineticLawName','X'))
sbioremovefromlibrary('kineticlaw','X')

X = sbioabstractkineticlaw("X","a*b*kf - c*kf/ke");
set(X,'SpeciesVariables',{'a' 'b' 'c'});
set(X,'ParameterVariables',{'kf' 'ke'});
sbioaddtolibrary(X)

addparameter(m,'kz',kz);
addparameter(m,'kf',kf);
kin = addkineticlaw(m.reaction(1),'X');
set(kin, 'parametervariablenames',{'kf' 'kz'});
setspecies(kin,'a','A');
setspecies(kin,'b','B');
setspecies(kin,'c','C');

conf = m.config;
set(conf.solveroptions, 'AbsoluteTolerance',1e-10);
set(conf,'stoptime',1e-1)

sim = sbiosimulate(m)
sbioplot(sim)

%% MassAction
clear all 

m = sbiomodel('model');

initial = 20

addspecies(m,'A',initial);
addspecies(m,'B',initial);
addspecies(m,'C',0);

addreaction(m,'A + B <-> C');

kf = 2e7;
kb = 2;

addparameter(m,'kb',kb);
addparameter(m,'kf',kf);
kin = addkineticlaw(m.reaction(1),'MassAction');
set(kin, 'parametervariablenames',{'kf' 'kb'});

conf = m.config;
set(conf.solveroptions, 'AbsoluteTolerance',1e-10);
runtime = 1e-8;
set(conf,'stoptime',runtime);

sim = sbiosimulate(m)

%%
A = sim.data(:,1);
e = 1;
index = find((A < (initial/2) + e) & (A > (initial/2) - e))
sbioplot(sim)
halftime = mean([sim.time(index(1)),sim.time(index(2))])
xline(halftime)
