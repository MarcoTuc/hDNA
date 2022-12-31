%% General Constants

%Avogadro Number 1/(mol)
NA      = 6.023e23;
%Nepero's constant 
gamma   = 0.57722;
%Boltzmann constant J/K
kboltz  = 1.38e-23;
%Gas constant (Kcal)/(K)(mol)
R_const = 1.987e-3;
%Planck constant (J)/(Hz)
kplanck = 6.626e-34;


%% Specific physical values

%Volume of the container
Volume = 1;

%Viscosity of water at 26°C (kg)/(m)(s) == (N)(s)/(m^2)
viscosity_h2o = 0.8701e-3;

%Diffusion Constant of a phospholipid molecule at 55% wetting in a membrane - Filippov2004 (dm^2)/(s)
phospholipid_diffusion = 9.32e-10;

%Radius of a cylinder representing a DNA double-helix 
dna_radius = 2e-9;

%steric factor of DNA renaturation hypothesized by christensen2001
steric_factor = 1/200;

%temperature 
temperature = 26;           % °C
T = 273.15 + temperature;   % °K

%% Formulas

%Stokes-Einstein solution diffusion coefficient
D_se = (kboltz*T)/(6*pi*viscosity_h2o*dna_radius);


%% Unit conversion

cal_joule = 4.1868;


%% Functions

function maxes = plotmax(data)
maxes = [];
for i = 1:length(data)
[max_N, indx_N] = max(data(i).Data(:,1));
tmax_N  = data(i).time(indx_N);
maxes = [maxes; max_N, tmax_N]; 
end
hold on
plot(maxes(:,2),maxes(:,1),'r+')
end
