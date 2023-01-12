clear all

kset = "explore";
modelname = "lida";
lida = buildmodel(modelname,kset);
v1 = lida.variant(1);
cursorvar = lida.variant(2);

species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];


%% Variant change of parameters

abstol = 1e-25;
reltol = 1e-13;
runtime_variant = 4e4;
timeunits = 'second';

f.setmodel_ode15s(lida,abstol,reltol,runtime_variant,timeunits);

v1.content = {{'parameter','kf_TR1','value',1.5e8}...
              {'parameter','kf_TL1','value',1.5e8}...
              {'parameter','kf_duplex','value',1.5e8}...
              {'parameter','kf_O1','value',1.5e8}};

data = sbiosimulate(lida,v1)
hold on
% sbioplot(data)
figure(1)
plot(data.time,data.data(:,12))
title('variant change')

%% Direct Change of parameters

kset = "explore";
modelname = "lida";
lida = buildmodel(modelname,kset);

abstol = 1e-25;
reltol = 1e-13;
runtime = 4e4
f.setmodel_ode15s(lida,abstol,reltol,runtime,timeunits);

data2 = sbiosimulate(lida)
hold on 
figure(2)
% sbioplot(sbio)
plot(data2.time,data2.data(:,12))
title('direct change')
