function [ e, a, b, c, thc, RHO, drho, PH, dph, TH, THi ] = calculateLensParameters( D, er, N )
%calculateLensParameters This function calculates parameters of untruncated
%elliptical lens
%   Detailed explanation goes here
%   Note: fix in-line documentation
    %% Calculate Eccintercity
    e = 1 / sqrt(er);
    %% Calculate Critical Angle
    thc = asin(e);
    %% Calculate Maximum Inclination Angle
    thmax = pi / 2 - thc;
    %% Calculate Minimum Radial Distance
    rmin = D / ( 2 * sin(thmax) );
    %% Calculate Lens Axes and Foci Distance
    a = rmin * ( 1 - e * cos(thmax) ) / (1 - e ^ 2);
    c = a * e;
    b = sqrt( a ^ 2 - c ^ 2);
    %% RHO and PHI of Cylindrical Coordinates
    rho = linspace(eps, D/2, N);
    drho = rho(2) - rho(1);
    ph = linspace(eps, 2 * pi, N);
    dph = ph(2) - ph(1);
    [ RHO, PH ] = meshgrid( rho, ph );
    %% Calculate Inclination Angle
    Z = a * sqrt( 1 - ( RHO / b ) .^ 2 ) + c;
    TH = atan( RHO ./ Z );
    %% Calculate Incident Angle
    THi = acos( ( 1 - e * cos(TH) ) ./ sqrt( 1 + e^2 - 2 * e * cos(TH) ) );
end