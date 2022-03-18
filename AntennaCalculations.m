close all;
clear;
clc;

%% Constants and Parameters
% Indecies
N = 200;                        % Number of sample points
% Field Parameters
f = 300 * 1e9;                  % Source frequency [Hz]
R = 1;                          % Radial distance [m]
p = [0 1 0];                    % Polarization of uniform current
J0 = 1;                         % Amplitude of uniform current
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

%% Theta and Phi-Components of Spherical Coordinates Outside Lense
thf = linspace(eps, pi / 2, N);
dthf = thf(2) - thf(1);
phf = linspace(eps, 2 * pi, N);
dphf = ph(2) - ph(1);
[ THf, PHf ] = meshgrid( thf, phf );

%% x, y, and z-Components of the wave-number
KX = k * sin(THf) .* cos(PHf);
KY = k * sin(THf) .* sin(PHf);
KZ = k * cos(THf);

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
[ Ef, Pfrad ] = calculateLensFeed( R, TH, PH, kd, Zd, n );
Ef = Ef * R / exp(-1j * kd * R);

%% Calculate Fresnel Transmission Coefficients
[ Tper, Tpar, THt ] = calculateFresnelTCoeff( THi, erd );

%% Calculate Equivalent Electric Current
[ J, M ] = calculateLensAperture( Ef, TH, PH, r, THt, THi, Tper, Tpar, ...
                                                                    Z, e );
plotCurrent( J, RHO, PH, 'J' );

%% Calculate Fourier Transform (FT) of Current
Jft = calculateCylFTCurrent( J, KX, KY, RHO, PH );
Jft = convertCylToCar(Jft, PHf);

%% Calculate Spectral Green's Function (SGF)
ej_SGF = calculateEJ_SGF( Z, k, KX, KY, KZ );

%% Calculate Electric Far-Field
E = calculateEFarfield( ej_SGF, Jft, k, R, THf, KZ );
E = convertCarToSph(E, THf, PHf);
plotFarfield(E, THf, PHf);
caxis([-40, 0]);
zlim([-150 0]);

%% Calculate Uniform Aperture Current Distribution Fourier Transform (FT)
Jun = circFTCurrent( k, J0, THf, D / 2, p );
Jun = convertSphToCar(Jun, THf, PHf);

%% Calculate Electric Far-Field of Uniform Aperture
Eun = calculateEFarfield( ej_SGF, Jun, k, R, THf, KZ );
Eun = convertCarToSph(Eun, THf, PHf);
Eun( isnan(Eun) ) = 0;
plotFarfield(Eun, THf, PHf);
caxis([-40, 0]);
zlim([-150 0]);

%% Plot Far-Field in 1D
% Define Theta from - Theta_max to Theta_max
thp = zeros( 1, 2 * size(THf, 2) );
thp( size(THf, 2) + 1 : end ) = THf(1, :);
thp( 1 : size(THf, 2) ) = - rot90( THf(1, :), 2 );
% Extract E field magnitude and uniform E field
Eth = zeros( 1, size(thp, 2) );
Eth( 1 : size(E, 2) ) = rot90( sqrt( abs( E(101, :, 1) ).^2 + ...
                    abs( E(101, :, 2) ).^2 + abs( E(101, :, 3) ).^2 ), 2 );
Eth( size(E, 2) + 1 : end ) = sqrt( abs( E(1, :, 1) ).^2 + ...
                             abs( E(1, :, 2) ).^2 + abs( E(1, :, 3) ).^2 );
Euth = zeros( 1, size(thp, 2) );
Euth( 1 : size(E, 2) ) = rot90( sqrt( abs( Eun(101, :, 1) ).^2 + ...
                 abs( Eun(101, :, 2) ).^2 + abs( Eun(101, :, 3) ).^2 ), 2);
Euth( size(E, 2) + 1 : end ) = sqrt( abs( Eun(1, :, 1) ).^2 + ...
                         abs( Eun(1, :, 2) ).^2 + abs( Eun(1, :, 3) ).^2 );
% Plot (normalized to uniform current field maximum)
figure();
plot(thp * 180 / pi, 20 * log10( Eth ) - max( 20 * log10( Euth ) ), ...
     'LineWidth', 3.0);
hold on;
plot(thp * 180 / pi, 20 * log10( Euth ) - max( 20 * log10( Euth ) ), ...
     '--', 'LineWidth', 3.0);
grid on;
xlabel('\theta [deg]');
ylabel('[dB]');
xlim([min(thp * 180 / pi) max(thp * 180 / pi)]);
xlim([-40 40]);
ylim([-40 0]);
legend('|E|', '|E_{un}|');
title('XZ Plane');
% Plot (normalized to own maximum)
figure();
plot(thp * 180 / pi, 20 * log10( Eth ) - max( 20 * log10( Eth ) ), ...
     'LineWidth', 3.0);
hold on;
plot(thp * 180 / pi, 20 * log10( Euth ) - max( 20 * log10( Euth ) ), ...
     '--', 'LineWidth', 3.0);
grid on;
xlabel('\theta [deg]');
ylabel('[dB]');
xlim([-40 40]);
ylim([-40 0]);
legend('|E|', '|E_{un}|');
title('XZ Plane');

%% Calculate Radiation Intensity and Radiation Power
Et = abs( E(:, :, 1) ) .^2 + abs( E(:, :, 2) ) .^2 + abs( E(:, :, 3) ) .^2;
U = Et * (R ^ 2) / (2 * Z);
Prad = sum( sum( U .* sin(THf) ) ) * dthf * dphf;

%% Calculate Directivity
D = 4 * pi * U / Prad;

%% Calculate Efficiency
Ae = Prad / Pfrad;

%% Calculate Gain
G = D * Ae;

%% Plot Antenna Parameters in 1D
% Extract directivity
Dth = zeros( 1, size(thp, 2) );
Dth( 1 : size(D, 2) ) = rot90( D(101, :), 2 );
Dth( size(G, 2) + 1 : end ) = D(1, :);
% Extract gain
Gth = zeros( 1, size(thp, 2) );
Gth( 1 : size(G, 2) ) = rot90( G(101, :), 2 );
Gth( size(G, 2) + 1 : end ) = G(1, :);
% Plot
figure();
plot(thp * 180 / pi, 10 * log10( Dth ), 'LineWidth', 3.0);
hold on;
plot(thp * 180 / pi, 10 * log10( Gth ), '--', 'LineWidth', 3.0);
grid on;
xlabel('\theta [deg]');
ylabel('[dB]');
xlim([min(thp * 180 / pi) max(thp * 180 / pi)]);
xlim([-40 40]);
ylim([-10 30]);
legend('D','G');
