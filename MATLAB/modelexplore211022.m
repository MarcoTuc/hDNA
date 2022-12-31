clear all 
buildmodel

%% Simulation Configuration

abstol = 1e-20;
reltol = 1e-13;
runtime = 2e6;
timeunits = 'second';

conf = lida.config;

set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);

set(conf, 'stoptime',runtime);
sbioaccelerate(lida);

%% Rates Exploration 1

exp1.name = 'kf_2d';
exp1.min  = 6;  %2eX
exp1.max  = 13; %2eX
exp1.num  = 40;

exp2.name = 'kf_3d';
exp2.min  = 6; %2eX
exp2.max  = 8; %2eX
exp2.num  = 40;

outputs = {'Nl'};
runtime = 2e10;
tic
data = FUNZ.explore_duplexequalized(lida, exp1, exp2, outputs, runtime);
toc
FUNZ.logallmaxsurf(data,exp1,exp2);


%% Single run Variant
i = 1
peaktime = zeros(1,2);
for kbd = [0.02, 2e2]
v1.content = {{'parameter','kf_lig','value',0.02}...
              {'parameter','kb_duplex','value',kbd}};
var = (sbiosimulate(lida,v1));
sbioplot(var);
peaktime(:,i) = FUNZ.peakplot(var,12);
FUNZ.variantlegend(v1)
i = i + 1; 
end

%% Surface Variant Exploration 1

v1.content = {{'parameter','kf_lig','value',0.002}...
              {'parameter','kb_duplex','value',0.03}...
              {'parameter','kb_oligo_bulk','value',2}};

exp1.name = 'kf_2d';
%logspaced
exp1.min  = 7;
exp1.max  = 12;
exp1.num  = 40;

exp2.name = 'kb_oligo_surf';
%logspaced
exp2.min  = 0;
exp2.max  = 5;
exp2.num  = 40;

outputs = {'Nl'};
runtime = 2e13;
tic
data = FUNZ.explore_duplexequalized(lida, v1, exp1, exp2, outputs, runtime);
toc

hold on
figure(1)
FUNZ.logallmaxsurf(data,exp1,exp2);
FUNZ.variantlegend(v1.content)

%% Surface Variant Exploration 2

v1.content = {{'parameter','kf_lig','value',0.02}...
              {'parameter','kb_duplex','value',3e1}...
              {'parameter','kb_oligo_bulk','value',2}};

exp1.name = 'kf_2d';
%logspaced
exp1.min  = 7;
exp1.max  = 12;
exp1.num  = 40;

exp2.name = 'kf_mix';
%logspaced
exp2.min  = 6;
exp2.max  = 8;
exp2.num  = 40;

outputs = {'Nl'};
runtime = 2e12;
tic
data = FUNZ.explore_duplexbulkequalized(lida, v1, exp1, exp2, outputs, runtime);
toc
hold on
figure(2)
FUNZ.logallmaxsurf(data,exp1,exp2);
FUNZ.variantlegend(v1.content)

%% Surface Variant Exploration 3


v1.content = {{'parameter','kf_2d','value',2e11}...
              {'parameter','kf_2d_duplex','value',2e11}...
              {'parameter','kf_3d','value',2e7}...
              {'parameter','kf_mix','value',2e8}};

exp1.name = 'kf_lig';
%logspaced
exp1.min  = -4;
exp1.max  = -1;
exp1.num  = 40;

exp2.name = 'kb_duplex';
%logspaced
exp2.min  = -3;
exp2.max  = 3;
exp2.num  = 40;

outputs = {'Nl'};
runtime = 2e12;
tic
data = FUNZ.explore_duplexbulkequalized(lida, v1, exp1, exp2, outputs, runtime);
toc

hold on
figure(3)
FUNZ.logallmaxsurf(data,exp1,exp2);
FUNZ.variantlegend(v1.content)



%% exp1
runtime = 2e7;
paramnames = {'kf_2d' 'kf_2d_duplex' 'kb_duplex'};
paramvalue = [2e7 2e7 2e2];
figure(1)
hold on
FUNZ.pointexplore(lida,paramnames,paramvalue,species_names,runtime);
legend(species_names);
% FUNZ.peakplot(point,12);
title('kf2dexplore');

%% exp2 - kset = 'mix-free'
paramnames = {'kf_2d' 'kf_2d_duplex' 'kf_3d' 'kf_mix' 'kb_oligo_surf' 'kb_oligo_bulk' 'kb_duplex' 'kf_lig'};
paramvalue = [1.5e12     1.5e12       2e7      2e7      7e5            1.5e2            3e-1        2e-2  ];
figure(2)
hold on
FUNZ.pointexplore(lida,paramnames,paramvalue,species_names,2e5);
legend(species_names);
% FUNZ.peakplot(point,12);
title('ratio 2d 3d explore');

%% exp2 - kset = 'mix-constrained'
paramnames = {'kf_2d'  'kf_3d' 'kf_mix' 'kb_oligo_surf' 'kb_oligo_bulk' 'kb_duplex' 'kf_lig'};
paramvalue = [1.5e12    2e7      2e7      7e5            1.5e2            3e-1        2e-2  ];
figure(3)
hold on
sbiosol = FUNZ.pointexplore(lida,paramnames,paramvalue,species_names,2e5);
legend(species_names);
% FUNZ.peakplot(point,12);
title('ratio 2d 3d explore');

%% cursor variant

cursorvar.content = {{'parameter','kf_2d','value',cursor2.Position(1)},...
                     {'parameter','kb_oligo_surf','value',cursor2.Position(2)},...
                     {'parameter','kf_lig','value',0.02}...
                     {'parameter','kb_duplex','value',0.03}...
                     {'parameter','kb_oligo_bulk','value',2}}

set(lida.config,'stoptime',2e3);
cursordata = sbiosimulate(lida,cursorvar)
sbioplot(cursordata)
cursorpeaktime = FUNZ.peakplot(cursordata,12)

%% att
for i = 1:15
sd_t = SD.Data(1:end-100,i);
sd_t1 = SD.Data(101:end,i);
figure(i)
plot(sd_t, sd_t1)
end

%% Conserved
SD = sbiosimulate(lida)

Tl = SD.Data(:,1);
R1 = SD.Data(:,2);
R2 = SD.Data(:,3);
TR1= SD.Data(:,4);
TR2= SD.Data(:,5);
Nl = SD.Data(:,6);
Tr = SD.Data(:,7);
L1 = SD.Data(:,8);
L2 = SD.Data(:,9);
TL1= SD.Data(:,10);
TL2= SD.Data(:,11);
Nr = SD.Data(:,12);
D  = SD.Data(:,13);
R  = SD.Data(:,14);
L  = SD.Data(:,15);

cons1 = R2 + TR2 + Nl + Tr + TL1 + TL2 + Nr + D + L;
cons2 = R1 + TR1 + Nl + Tr + TL1 + TL2 + Nr + D + L;
cons3 = Tl + TR1 + TR2 + Nl + L2 + TL2 + Nr + D + L;
cons4 = Tl + TR1 + TR2 + Nl + L1 + TL1 + Nr + D + R;
figure(1)
plot(SD.Time,cons1)
ylim([0 2e-4])
figure(2)
plot(SD.Time,cons2)
ylim([0 2e-4])
figure(3)
plot(SD.Time,cons3)
ylim([0 2e-4])
figure(4)
plot(SD.Time,cons4)
ylim([0 2e-4])


rcons1 = - Tl + R2 - TR1 + Tr - L2 + TL1;
rcons2 = - Tl + R1 - TR2 + Tr - L1 + TL2;
rcons3 = 2*Tl - R2 + 2*TR1 + TR2 + Nl - Tr + L1 + L2 + Nr + D + R;
rcons4 = 2*Tl - R1 + TR1 + 2*TR2 + Nl - Tr + L1 + L2 + Nr + D + R;
figure(5)
plot(SD.Time,rcons1)
ylim([0 2e-4])
figure(6)
plot(SD.Time,rcons2)
ylim([0 2e-4])
figure(7)
plot(SD.Time,rcons3)
ylim([0 2e-4])
figure(8)
plot(SD.Time,rcons4)
ylim([0 2e-4])
