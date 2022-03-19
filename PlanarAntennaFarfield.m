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

%% Calculate Electric Far-Field Inside Lens
[ Ef, Pfrad ] = calculateLensFeed( R, TH, PH, kd, Zd, n );
plotFarfield( Ef, TH, PH );
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
xticks((-0.5 : 0.1 : 0.5));
yticks((-0.5 : 0.1 : 0.5));

%% Plot Far-Field in 1D
% Define Theta from - Theta_max to Theta_max
thp = zeros( 1, 2 * size(TH, 2) );
thp( size(TH, 2) + 1 : end ) = TH(1, :);
thp( 1 : size(TH, 2) ) = - rot90( TH(1, :), 2 );
% Extract E field magnitude
Eth = zeros( 1, size(thp, 2) );
Eth( 1 : size(Ef, 2) ) = rot90( sqrt( abs( Ef(251, :, 1) ).^2 + ...
                  abs( Ef(251, :, 2) ).^2 + abs( Ef(251, :, 3) ).^2 ), 2 );
Eth( size(Ef, 2) + 1 : end ) = sqrt( abs( Ef(1, :, 1) ).^2 + ...
                           abs( Ef(1, :, 2) ).^2 + abs( Ef(1, :, 3) ).^2 );
% Plot (normalized to own maximum)
figure();
plot(thp * 180 / pi, 20 * log10( Eth ) - max( 20 * log10( Eth ) ), ...
     'LineWidth', 3.0);
grid on;
xlabel('\theta [deg]');
ylabel('|E| [dB]');
xlim([min(thp * 180 / pi) max(thp * 180 / pi)]);
ylim([-40 0]);
title('XZ Plane');
