close all;
clear;
clc;

%% Constants and Parameters
% Indecies
N = 500;                        % Number of sample points
% Field Parameters
f = 300 * 1e9;                  % Source frequency [Hz]
R = 1;                          % Radial distance [m]
n = 3;
% Medium
er = 1;                         % Relative permittivity
erd = 11.9;                     % Relative permittivity inside lens
c = physconst('LightSpeed');    % Speed of light [m/s]
e0 = 8.8541878128 * 1e-12;      % Permittivity of free space [F/m]
u0 = 4 * pi * 1e-7;             % Permeability of free space [H/m]

%% Parameters
wlen = c / f;                   % Wavelength [m]
k = 2*pi / wlen;                % Magnitude of wave number [rad/m]
kd = k * sqrt(erd);             % Magnitude of wave number inside lens [rad/m]
Z = sqrt( u0 / (e0 * er) );     % Wave impedance [Ohm]
Zd = Z / sqrt(erd);             % Wave impedance inside lens [Ohm]

%% Calculate Crticial Angle
thc = asin( 1 / sqrt(erd) );

%% Theta and Phi-Components of Spherical Cooridnates
th = linspace(eps, thc, N);
dth = th(2) - th(1);
ph = linspace(eps, 2 * pi, N);
dph = ph(2) - ph(1);
[ TH, PH ] = meshgrid(th, ph);

%% Calculate Fresnel Transmission Coefficients
[ Tper, Tpar, THt ] = calculateFresnelTCoeff( TH, erd );

%% Calculate Transmitted Power Ratio
[ Prper, Prpar ] = calculateTPowerRatio( Tper, Tpar, TH, THt, erd );

%% Plot
figure();
plot( th * 180 / pi, Prper(1, :), 'LineWidth', 3.0 );
hold on;
plot( th * 180 / pi, Prpar(1, :), 'LineWidth', 3.0 );
grid on;
xlabel('\theta_{i} [deg]'); 
legend('P_{t}^{TE} / P_{i}^{TE}', 'P_{t}^{TM} / P_{t}^{TM}');
xlim([0 max(th * 180 / pi)]);
