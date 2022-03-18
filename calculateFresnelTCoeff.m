function [ Tper, Tpar, THt ] = calculateFresnelTCoeff( THi, er )
%calculateFresnelTCoeff This function calculates the Fresnel transmission
%coefficients and the transmission angle
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Define Wave Impedance
    Z0 = 376.73031346177;   % Wave impedance in free space [Ohm]
    Z = Z0 / sqrt(er);      % Wave impedance in lens medium [Ohm]
    %% Calculate Transmission Angle
    THt = asin( sqrt(er) * sin(THi) );
    %% Calculate Transmission Coefficients
    Tper = 2 * Z0 * cos(THi) ./ ( Z0 * cos(THi) + Z * cos(THt) );
    Tpar = 2 * Z0 * cos(THi) ./ ( Z0 * cos(THt) + Z * cos(THi) );
end