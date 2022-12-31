clear all

constants2D

%% Initial conditions/concentrations

y0_oligo = 1e-4;
y0_temp  = 1e-9;

initial = [];

initial(1)  = y0_oligo; %Or1
initial(2)  = y0_oligo; %Or2
initial(3)  = y0_oligo; %Ol1
initial(4)  = y0_oligo; %Ol2   
initial(5) = 0;         %orl_1
initial(6) = 0;         %orl_2
initial(7) = y0_temp;   %Tl - The Guy    
initial(8) = 0;         %Tr
initial(9) = 0;         %Td         
initial(10) = 0;        %Tl_Or1
initial(11) = 0;        %Tl_Or2
initial(12) = 0;        %Tl_Or1_Or2
initial(13) = 0;        %Tr_Ol1
initial(14) = 0;        %Tr_Ol2
initial(15) = 0;        %Tr_Ol1_Ol2

%% Mass Action ODEs

LIDA = @(t,y) [
%Or1 y(1)
kb_or1_tleft*(y(10)+y(12)) - kf_or1_tleft*(y(1))*(y(7)+y(11)) - kf_orl_1*(y(1))*(y(3)) + kb_orl_1*y(5);
%Or2 y(2)
kb_or2_tleft*(y(11)+y(12)) - kf_or2_tleft*(y(2))*(y(7)+y(10)) - kf_orl_2*(y(2))*(y(4)) + kb_orl_2*(y(6));
%Ol1 y(3)
kb_ol1_tright*(y(13)+y(15)) - kf_ol1_tright*(y(3))*((y(8)+y(14))) - kf_orl_1*(y(1))*(y(3)) + kb_orl_1*y(5);
%Ol2 y(4)
kb_ol2_tright*(y(14)+y(15)) - kf_ol2_tright*(y(4))*((y(8)+y(13))) - kf_orl_2*(y(2)*y(4)) + kb_orl_2*(y(6));
%orl_1 y(5)
kf_orl_1*(y(1)*y(3)) - kb_orl_1*y(5);
%orl_2 y(6)
kf_orl_2*(y(2)*y(4)) - kb_orl_2*y(6);
%Tl y(7)
kb_duplex*(y(9)) + kb_or2_tleft*(y(11)) + kb_or1_tleft*(y(10)) - kf_duplex*(y(8))*(y(7)) - kf_or1_tleft*(y(1))*(y(7)) - kf_or2_tleft*(y(2))*(y(7));
%Tr y(8)
kb_duplex*(y(9)) + kb_ol1_tright*(y(13)) + kb_ol2_tright*(y(14)) - kf_duplex*(y(7))*(y(8)) - kf_ol1_tright*(y(3))*(y(8)) - kf_ol2_tright*(y(4))*(y(8));
%Td y(9)
kf_duplex*(y(7))*(y(8)) + kf_repair*(y(12)+y(15)) - kb_duplex*(y(9));
%Tl_or1 y(10)
kf_or1_tleft*(y(7))*(y(1)) + kb_or2_tleft*(y(12)) - kf_or2_tleft*(y(10))*(y(2)) - kb_or1_tleft*(y(10));
%Tl_or2 y(11)
kf_or2_tleft*(y(7))*(y(2)) + kb_or1_tleft*(y(12)) - kf_or1_tleft*(y(11))*(y(1)) - kb_or2_tleft*(y(11));
%Tl_or1_or2 y(12)
kf_or1_tleft*(y(11))*(y(1)) + kf_or2_tleft*(y(10))*(y(2)) - kf_repair*(y(12)) - kb_or1_tleft*(y(12)) - kb_or2_tleft*(y(12));
%Tr_ol1 Y(13)
kf_ol1_tright*(y(8))*(y(3)) + kb_ol2_tright*(y(15)) - kf_ol2_tright*(y(13))*(y(4)) - kb_ol1_tright*(y(13));
%Tr_ol2 y(14)
kf_ol2_tright*(y(8))*(y(4)) + kb_ol1_tright*(y(15)) - kf_ol1_tright*(y(14))*(y(3)) - kb_ol2_tright*(y(14));
%Tr_ol1_ol2 y(15)
kf_ol1_tright*(y(14))*(y(3)) + kf_ol2_tright*(y(13))+(y(4)) - kf_repair*(y(15)) - kb_ol2_tright*(y(15)) - kb_ol1_tright*(y(15));
    ];

options = odeset('RelTol',1e-13,'AbsTol',1e-30);
runtime = 2e6;

tic
[t,c] = ode23s(LIDA,(0:runtime),initial);
toc

figure(1)
hold on
plot(t/60,c)

legend('O_{r,1}','O_{r,2}','O_{l,1}','O_{l,2}','O_{rl,1}','O_{rl,2}', ...
    'T_l','T_r','T_d', ...
    'T_lO_{r,1}','T_lO_{r,2}','T_lO_{r,1}O_{r,2}', ...
    'T_rO_{l,1}','T_rO_{l,2}','T_rO_{l,1}O_{l,2}')
