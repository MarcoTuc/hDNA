kset = "mix-constrained";
modelname = "lida_constr";
lida_constr = buildmodel(modelname,kset) 

species_names   = ["Tl", "R1", "R2", "TR1", "TR2", "Nl", "Tr", "L1", "L2", "TL1", "TL2", "Nr", "D", "R", "L"];

v1 = lida_constr.variant(1)
cursorvar = lida_constr.variant(2)

%% Simulation Configuration 

abstol = 1e-20;
reltol = 1e-13;
runtime = 2e6;
timeunits = 'second';

f.setmodel_ode15s(lida_constr,abstol,reltol,runtime,timeunits);

%% Conserved Moieties 
SD = sbiosimulate(lida_constr)
plot = false;
% Assigning meaningful names to variables
Tl = SD.Data(:,1); R1 = SD.Data(:,2); R2 = SD.Data(:,3); TR1= SD.Data(:,4); TR2= SD.Data(:,5); Nl = SD.Data(:,6);
Tr = SD.Data(:,7); L1 = SD.Data(:,8); L2 = SD.Data(:,9); TL1= SD.Data(:,10); TL2= SD.Data(:,11); Nr = SD.Data(:,12);
D  = SD.Data(:,13); R  = SD.Data(:,14); L  = SD.Data(:,15);

% Conservation laws obtained through sbioconsmoiety(lida_constr,'semipos','p')
conserved_semipos = sbioconsmoiety(lida_constr,'semipos','p');
cons1 = R2 + TR2 + Nl + Tr + TL1 + TL2 + Nr + D + L;
cons2 = R1 + TR1 + Nl + Tr + TL1 + TL2 + Nr + D + R;
cons3 = Tl + TR1 + TR2 + Nl + L2 + TL2 + Nr + D + L;
cons4 = Tl + TR1 + TR2 + Nl + L1 + TL1 + Nr + D + R;

if plot == true
figure(1); plot(SD.Time,cons1); ylim([0 2e-4]); title("cons1",conserved_semipos(1))
figure(2); plot(SD.Time,cons2); ylim([0 2e-4]); title("cons2",conserved_semipos(2))
figure(3); plot(SD.Time,cons3); ylim([0 2e-4]); title("cons3",conserved_semipos(3))
figure(4); plot(SD.Time,cons4); ylim([0 2e-4]); title("cons4",conserved_semipos(4))
end
% Conservation laws obtained thorugh sbioconsmoiety(lida_constr,'rreduce','p')
conserved_rreduce = sbioconsmoiety(lida_constr,'rreduce','p');
rcons1 = - Tl + R2 - TR1 + Tr - L2 + TL1;
rcons2 = - Tl + R1 - TR2 + Tr - L1 + TL2;
rcons3 = 2*Tl - R2 + 2*TR1 + TR2 + Nl - Tr + L1 + L2 + Nr + D + R; 
rcons4 = 2*Tl - R1 + TR1 + 2*TR2 + Nl - Tr + L1 + L2 + Nr + D + L; 
if plot == true
figure(5); plot(SD.Time,rcons1); ylim([0 2e-4]); title("rcons1",conserved_rreduce(1))
figure(6); plot(SD.Time,rcons2); ylim([0 2e-4]); title("rcons2",conserved_rreduce(2))
figure(7); plot(SD.Time,rcons3); ylim([0 2e-4]); title("rcons3",conserved_rreduce(3))
figure(8); plot(SD.Time,rcons4); ylim([0 2e-4]); title("rcons4",conserved_rreduce(4))
end

%% Extended kb_duplex vs kf_lig 

fwd_2d = 4e8;
fwd_3d = 2e7;

v1.content = {{'parameter','kf_2d','value',fwd_2d}...
              {'parameter','kf_3d','value',fwd_3d}...
              {'parameter','kf_mix','value',fwd_3d}};

exp1.name = 'kf_lig';
%logspaced
exp1.min  = -5;
exp1.max  = 0;
exp1.num  = 40;

exp2.name = 'kb_duplex';
%logspaced
exp2.min  = -10;
exp2.max  = 8;
exp2.num  = 40;

outputs = species_names;
runtime = 2e9;
tic
data = f.explore(lida_constr, v1, exp1, exp2, outputs, runtime);
toc

hold on
figure(1)
f.logallmaxsurf_mod(data,12,exp1,exp2);
legendcontent = {{'parameter','kf_2d','value',fwd_2d}...
                {'parameter','kf_3d = kf_mix','value',fwd_3d}};
f.variantlegend(legendcontent)


%% cursor variant green
cursor1 = tip_green;
cursorvar.content = {   {'parameter','kf_lig','value',cursor1.Position(1)},...
                        {'parameter','kb_duplex','value',cursor1.Position(2)}...
                        {'parameter','kf_2d','value',fwd_2d}...
                        {'parameter','kf_3d','value',fwd_3d}...
                        {'parameter','kf_mix','value',fwd_3d}};
cursor_runtime = 2e13;
set(lida_constr.config,'stoptime',cursor_runtime);
cursordata = sbiosimulate(lida_constr,cursorvar);
sbioplot(cursordata)
title('green')
cursorpeaktime_green = f.peakplot(cursordata,12)

%% cursor variant deep blue
cursor2 = tip_blue;
cursorvar.content = {   {'parameter','kf_lig','value',cursor2.Position(1)},...
                        {'parameter','kb_duplex','value',cursor2.Position(2)}...
                        {'parameter','kf_2d','value',fwd_2d}...
                        {'parameter','kf_3d','value',fwd_3d}...
                        {'parameter','kf_mix','value',fwd_3d}}
cursor_runtime = 2e11;
set(lida_constr.config,'stoptime',cursor_runtime);
cursordata = sbiosimulate(lida_constr,cursorvar);
sbioplot(cursordata)
title('deepblue')
cursorpeaktime_blue = f.peakplot(cursordata,12)



%% Extended kf_2d searching for product inhibition

fwd_3d = 2e7;

v1.content = {{'parameter','kf_lig','value',2e-2}...
              {'parameter','kf_3d','value',fwd_3d}...
              {'parameter','kf_mix','value',fwd_3d}};

exp1.name = 'kf_2d';
%logspaced
exp1.min  = 6;
exp1.max  = 14;
exp1.num  = 20;

exp2.name = 'kb_duplex';
%logspaced
exp2.min  = -10;
exp2.max  = 8;
exp2.num  = 20;

outputs = species_names;
runtime = 2e10;
tic
data = f.explore(lida_constr, v1, exp1, exp2, outputs, runtime);
toc

hold on
figure(1)
f.maxsurf_mod(data,12,exp1,exp2,true);
legendcontent = {{'parameter','kf_lig','value',2e-2}...
              {'parameter','kf_3d','value',fwd_3d}}
f.variantlegend(legendcontent,8)

%% Cursorcheck

for i = 1:length(cursor_data)
cursor = cursor_data(i);
variantcontent = {{'parameter','kf_2d','value',cursor.Position(1)}...
                     {'parameter','kb_duplex','value',cursor.Position(2)}...
                     {'parameter','kf_lig','value',2e-2}...
                     {'parameter','kf_3d','value',fwd_3d}...
                     {'parameter','kf_mix','value',fwd_3d}};

cursorvar.content = variantcontent;
cursor_runtime = 3*cursor.Position(3);
set(lida_constr.config,'stoptime',cursor_runtime);
cursordata = sbiosimulate(lida_constr,cursorvar);
title_string = join(string(['cursor',' ',string(i)]));
sbioplot(cursordata);
title(title_string);
cursorpeaktime = f.peakplot(cursordata,12);
check = join({sprintf('%.3e',cursorpeaktime), sprintf('%.3e',cursor.Position(3))},newline);
annotation('textbox',[.75, .2, 0, 0],'String',check,'Interpreter','None','FitBoxToText','on','Edgecolor','None');
f.variantlegend(variantcontent(1:4),8);

end