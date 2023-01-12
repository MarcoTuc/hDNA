bimo = sbiomodel('bimo')

y0_R = 1
y0_O = 0

addspecies(bimo,'R',y0_R)
addspecies(bimo,'O',y0_O)

addreaction(bimo,'R + R <-> O')

kf = 2e7
kb = 2e7

addkineticlaw(bimo.reaction(1),'MassAction')
addparameter(bimo.reaction(1).kineticlaw,'Kf',kf)
addparameter(bimo.reaction(1).kineticlaw,'Kb',kb)
set(bimo.reaction(1).kineticlaw,'parametervariablenames',{'Kf' 'Kb'})

conf = bimo.config;

set(conf, 'solvertype','ode15s');
set(conf, 'stoptime',4e-7);
set(conf.solveroptions, 'absolutetolerance',1e-13);
set(conf.solveroptions, 'relativetolerance',1e-10);

SD = sbiosimulate(bimo)

sbioplot(SD)