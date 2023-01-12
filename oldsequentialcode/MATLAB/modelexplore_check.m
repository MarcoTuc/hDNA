kset = "explore";
modelname = "lida_check";
lida_check = buildmodel(modelname,kset)

species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];

%% Simulation Configuration 

abstol = 1e-20;
reltol = 1e-13;
runtime = 2e6;
timeunits = 'second';

conf = lida_check.config;

set(conf, 'solvertype','ode15s');
set(conf.solveroptions, 'absolutetolerance',abstol);
set(conf.solveroptions, 'relativetolerance',reltol);
set(conf, 'timeunits',timeunits);

set(conf, 'stoptime',runtime);
sbioaccelerate(lida_check);



%% Extended kb_duplex vs kf_lig 

fwd_2d = 2e7;
fwd_3d = 2e7;
kb_oligo_surf = 6e5;
kb_oligo_bulk = 120;
kb_duplex = 0.0072;
kf_lig = 0.02;

v1.content = {  {'parameter','kf_TR1','value',fwd_2d}...
              	{'parameter','kf_TR2','value',fwd_3d}...
                {'parameter','kf_TL1','value',fwd_2d}...
                {'parameter','kf_TL2','value',fwd_3d}...
                {'parameter','kb_TR1','value',kb_oligo_surf}...
              	{'parameter','kb_TR2','value',kb_oligo_bulk}...
                {'parameter','kb_TL1','value',kb_oligo_surf}...
                {'parameter','kb_TL2','value',kb_oligo_bulk}...
                {'parameter','kf_lig','value',kf_lig}...
                {'parameter','kf_duplex','value',fwd_2d}...
                {'parameter','kb_duplex','value',kb_duplex}...
                {'parameter','kf_O1','value',fwd_2d}...
                {'parameter','kf_O2','value',fwd_3d}};

exp1.name = 'kf_2d';
exp1.min  = 6; %logspaced
exp1.max  = 14;
exp1.num  = 20;

exp2.name = 'kb_duplex';
exp2.min  = -10; %logspaced
exp2.max  = 8;
exp2.num  = 20;

outputs = species_names;
runtime = 2e10;
tic
data = f.explore_check(lida_check, v1, exp1, exp2, outputs, runtime);
toc

hold on
figure(1)
f.maxsurf_mod(data,12,exp1,exp2,true);
legendcontent = {{'parameter','kf_lig','value',2e-2}...
              {'parameter','kf_3d','value',fwd_3d}}
f.variantlegend(legendcontent,8)

