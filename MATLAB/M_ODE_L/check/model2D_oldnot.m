

%% Initial conditions/concentrations
initial = [];

initial(1)  = 1e-4;         %O1T+
initial(2)  = 1e-4;         %O2T+ 
initial(3)  = 1e-4;         %O1T-
initial(4)  = 1e-4;         %O2T-
initial(5) = 0;             % tminuso1tplus
initial(6) = 0;             % tminuso2tplus
initial(7) = 0;             % tminuso1tpluso2tplus\
initial(8) = 0;             % tplustminus
initial(9) = 0;             % tpluso1tminus
initial(10) = 0;            % tpluso2tminus
initial(11) = 0;            % tpluso1tminuso2tminus
initial(12) = 1e-9;         % tminus
initial(13) = 0;            % tplus
initial(14) = 0;            % o1pluso1minus
initial(15) = 0;            % o2pluso2minus


%% Specification of the System
deq1=@(t,y) [

    kminuso1tplus*(y(5)+y(7))-kpluso1tplus*(y(1))*(y(12)+y(6))-kpluso1tpm*(y(1)*y(3))+kminuso1tpm*(y(14)); %[1]

    kminuso2tplus*(y(6)+y(7))-kpluso2tplus*(y(2))*(y(12)+y(5))-kpluso2tpm*(y(2)*y(4))+kminuso2tpm*(y(15)); %[2]

    kminuso1tminus*(y(9)+y(11))-kpluso1tminus*(y(13)*y(3)+y(10)*y(3))-kpluso1tpm*(y(1)*y(3))+kminuso1tpm*y(14); %[3]

    kminuso2tminus*(y(10)+y(11))-kpluso2tminus*(y(13)*y(4)+y(9)*y(4))-kpluso2tpm*(y(2)*y(4))+kminuso2tpm*y(15); %[4]

    kpluso1tplus*(y(12)*y(1))+kminuso2tplus*y(7)-kminuso1tplus*y(5)-kpluso2tplus*(y(5)*y(2)); %[5]

    kpluso2tplus*(y(12)*y(2))+kminuso1tplus*y(7)-kminuso2tplus*y(6)-kpluso1tplus*(y(6)*y(1)); %[6]

    kpluso1tplus*(y(6)*y(1))+kpluso2tplus*(y(5)*y(2))-kminuso1tplus*y(7)-kminuso2tplus*y(7)-klplus*y(7); %[7]

    klplus*y(7)+klminus*y(11)+ktplus*(y(12)*y(13))-ktminus*y(8); %[8]

    kpluso1tminus*(y(13)*y(3))+kminuso2tminus*y(11)-kminuso1tminus*y(9)-kpluso2tminus*(y(9)*y(4)); %[9]

    kpluso2tminus*(y(13)*y(4))+kminuso1tminus*y(11)-kminuso2tminus*y(10)-kpluso1tminus*(y(10)*y(3)); %[10]

    kpluso1tminus*(y(10)*y(3))+kpluso2tminus*(y(9)*y(4))-kminuso1tminus*y(11)-kminuso2tminus*y(11)-klminus*y(11); %[11]

    kminuso1tplus*y(5)+kminuso2tplus*y(6)+ktminus*y(8)-kpluso1tplus*(y(12)*y(1))-kpluso2tplus*(y(12)*y(2))-ktplus*(y(12)*y(13)); %[12]

    kminuso1tminus*y(9)+kminuso2tminus*y(10)+ktminus*y(8)-kpluso1tminus*(y(13)*y(3))-kpluso2tminus*(y(13)*y(4))-ktplus*(y(12)*y(13)); %[13]

    kpluso1tpm*(y(1)*y(3)) - kminuso1tpm*y(14); %[14]

    kpluso2tpm*(y(2)*y(4)) - kminuso2tpm*y(15); %[15]

    ];

%% Numerical Integration
options= odeset('RelTol',1e-9,'AbsTol',1e-12,'Refine',100);
runtime = 4e4;
step = runtime*10^-1;
tic
[t,sol] = ode15s(deq1,(0:runtime),initial,options);
toc

%% plot
figure(3)
hold on
plot(t,sol(:,11))
% legend('O_1^{T+}','O_2^{T+}','O_1^{T-}','O_2^{T-}','T^{-}O_1^{T+}', ...
%     'T^{-}O_2^{T+}','T^{-}O_1^{T+}O_2^{T+}','T^{-}T^{+}', ...
%     'T^{+}O_1^{T-}','T^{+}O_2^{T-}','T^{+}O_1^{T-}O_2^{T-}', ...
%     'T^{-}','T^{+}', ...
%     'O_1^{T+}O_1^{T-}','O_2^{T+}O_2^{T-}')
title('ODE model')


