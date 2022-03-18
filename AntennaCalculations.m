close all;
clear;
clc;

%% Constants and Parameters
% Indecies
N = 500;                        % Number of sample points
% Field Parameters
f = 300 * 1e9;                  % Source frequency [Hz]
R = 1;                          % Radial distance [m]
% Lens
n = 3;
th0 = 50 * pi / 180;            % Maximum inclination angle in lense [rad]
erd = 11.9;                     % Relative permittivity of lens medium
% Medium
er = 1;                         % Relative permittivity
c = physconst('LightSpeed');    % Speed of light [m/s]
e0 = 8.8541878128 * 1e-12;      % Permittivity of free space [F/m]
u0 = 4 * pi * 1e-7;             % Permeability of free space [H/m]

%% Parameters
% Wave
wlen = c / f;                   % Wavelength [m]
k = 2*pi / wlen;                % Magnitude of wave number [rad/m]
kd = k * sqrt(erd);             % Magnitude of wave number inside lens [rad/m]
% Wave impedance
Z = sqrt( u0 / (e0 * er) );     % Wave impedance [Ohm]
Zd = Z / sqrt(erd);             % Wave impedance inside lens [Ohm]
% Lense
D = 10 * wlen;                  % Diameter of lens [m]
e = 1 / sqrt(erd);              % Eccintercity

%% RHO and PHI of Cylindrical Coordinates
rho = linspace(eps, D/2, N);
drho = rho(2) - rho(1);
ph = linspace(eps, 2 * pi, N);
dph = ph(2) - ph(1);
[ RHO, PH ] = meshgrid( rho, ph );

%% Calculate Minimum Radial Distance
rmin = D / ( 2 * sin(th0) );

%% Calculate Lens Axes and Foci Distance
a = rmin * ( 1 - e * cos(th0) ) / (1 - e ^ 2);
c = a * e;
b = sqrt( a ^ 2 - c ^ 2);

%% Z of Cylindrical Coordinates
z = a * sqrt( 1 - ( RHO / b ) .^ 2 ) + c;

%% Calculate Inclination Angle
TH = atan( RHO ./ z );

%% Calculate Incident Angle
THi = acos( ( 1 - e * cos(TH) ) ./ sqrt( 1 + e^2 - 2 * e * cos(TH) ) );

%% Calculate Radial Distance
r = a * ( 1 - e^2 ) ./ ( 1 - e * cos(TH) );

%% Calculate Feed Electric Far-Field of Feed
[ Ef, ~ ] = calculateLensFeed( R, TH, PH, kd, Zd, n );
Ef = Ef * R / exp(-1j * kd * R);

%% Calculate Fresnel Transmission Coefficients
[ Tper, Tpar, THt ] = calculateFresnelTCoeff( THi, erd );

%% Calculate Equivalent Electric Current
[ J, M ] = calculateLensAperture( Ef, TH, PH, r, THt, THi, Tper, Tpar, ...
                                                                    Z, e );
plotCurrent( J, RHO, PH, 'J' );
