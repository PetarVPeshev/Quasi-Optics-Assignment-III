function [ E, Prad ] = calculateLensFeed( R, TH, PH, k, Z, n )
%calculateLensFeed This function calculates the electric far-field of a
%planar lens feed antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Extract Element Step
    dth = TH(1, 2) - TH(1, 1);
    dph = PH(2, 1) - PH(1, 1);
    %% Calculate Electric Far-Field
    E = zeros([ size(TH), 3 ]);
    E(:, :, 2) = - ( cos(TH) .^ n ) .* sin(PH) * exp(-1j * k * R) / R;
    E(:, :, 3) = ( cos(TH) .^ n ) .* cos(PH) * exp(-1j * k * R) / R;
    %% Calculate Total Electric Far-Field
    Et = abs( E(:, :, 1) ) .^2 + abs( E(:, :, 2) ) .^2 + ...
         abs( E(:, :, 3) ) .^2;
    %% Calculate Radiation Intensity
    U = Et * (R ^ 2) / (2 * Z);
    %% Calculate Radiated Power
    Prad = sum( sum( U .* sin(TH) ) ) * dth * dph;
end