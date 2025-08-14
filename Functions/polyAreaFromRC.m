function A = polyAreaFromRC(l, rc)
%==========================================================================
% polyAreaFromRC: Computes the area of a regular polygon in the plane
%                 from its circumradius.
%
% Inputs:
%   l   - Number of sides (integer, l >= 3).
%   rc  - Circumradius (same units as the desired linear measure).
%
% Output:
%   A   - Polygon area (square of the input length units).
%
% Formula Reference:
%   For a regular l-gon with circumradius rc:
%   A = (l/2) * rc^2 * sin(2*pi/l).
%==========================================================================

    % Compute area from the circumradius using the closed-form expression
    A = (l / 2) * rc^2 * sin(2 * pi / l);

end
