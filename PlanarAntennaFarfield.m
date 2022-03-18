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

%% Theta and Phi-Components of Spherical Cooridnates
th = linspace(eps, pi / 2, N);
dth = th(2) - th(1);
ph = linspace(eps, 2 * pi, N);
dph = ph(2) - ph(1);
[ TH, PH ] = meshgrid(th, ph);

%% x, y, and z-Components of Wave Number for Feed
KX = k * sin(TH) .* cos(PH);
KY = k * sin(TH) .* sin(PH);
kz = -1j * sqrt( -(k^2 - KX.^2 - KY.^2) );

%% Calculate Electric Far-Field Inside Lens
[ Ef, Pfrad ] = calculateLensFeed( R, TH, PH, kd, Zd, n );
plotFarfield( Ef, TH, PH );
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
xticks((-0.5 : 0.1 : 0.5));
yticks((-0.5 : 0.1 : 0.5));
