function [ J, M ] = calculateLensAperture( E, TH, PH, R, THt, THi, Tper, Tpar, Z, e )
%calculateLensAperture This function calculates the equivalent aperutre
%electric and magnetic current distributions of lens antenna
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Spreading Factor
    S = sqrt( ( cos(THt) .* ( e * cos(TH) - 1 ) ) ./ ...
              ( cos(THi) .* ( e - cos(TH) ) ) );
    %% Calculate Equivalent Electric Current Distribution
    J = zeros([ size(PH), 3 ]);
    J(:, :, 1) = - 2 * ( Tpar .* E(:, :, 3) .* cos(PH) - ...
                         Tper .* E(:, :, 2) .* sin(PH) ) .* S ./ (Z * R);
    J(:, :, 2) = - 2 * ( Tpar .* E(:, :, 3) .* sin(PH) + ...
                         Tper .* E(:, :, 2) .* cos(PH) ) .* S ./ (Z * R);
    %% Calculate Equivalent Magnetic Current Distirbution
    M = zeros([ size(PH), 3 ]);
    M(:, :, 1) = 2 * ( Tpar .* E(:, :, 3) .* cos(PH) - ...
                       Tper .* E(:, :, 2) .* sin(PH) ) .* S ./ R;
    M(:, :, 2) = 2 * ( Tpar .* E(:, :, 3) .* sin(PH) + ...
                       Tper .* E(:, :, 2) .* cos(PH) ) .* S ./ R;
end