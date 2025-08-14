function ga = goldenAngle()
%==========================================================================
% goldenAngle: Computes the golden angle in radians.
%
% Output:
%   ga   - Golden angle (radians).
%
% Definition:
%   Based on the golden ratio φ = (1 + sqrt(5)) / 2, the golden angle is:
%   ga = 2*pi*(1 - 1/φ) = pi*(3 - sqrt(5)).
%==========================================================================

    % Compute the golden angle in radians
    ga = pi * (3 - sqrt(5));  % From the golden ratio relationship

end
