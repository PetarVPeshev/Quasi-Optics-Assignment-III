function [ Prper, Prpar ] = calculateTPowerRatio( Tper, Tpar, THi, THt, er )
%calculateTPowerRatio This function calculates the transmitted power ratio
%of perpendicular and parallel polarizations
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Define Wave Impedance
    Z0 = 376.73031346177;   % Wave impedance in free space [Ohm]
    Z = Z0 / sqrt(er);      % Wave impedance in lens medium [Ohm]
    %% Calculate Transmitted Power Ratio
    Prper = ( abs( Tper ) .^ 2 ) .* Z .* cos(THt) ./ ( Z0 * cos(THi) );
    Prpar = ( abs( Tpar ) .^ 2 ) .* Z .* cos(THt) ./ ( Z0 * cos(THi) );
end