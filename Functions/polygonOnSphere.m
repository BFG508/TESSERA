function [P] = polygonOnSphere(center, rc, l, phase, R)
%==========================================================================
% polygonOnSphere: Builds a regular l-gon on the surface of a sphere of
%                  radius R, centered around the direction given by 'center'.
%
% Inputs:
%   center  - 3-element vector pointing to the desired polygon center
%             direction (will be normalized internally).
%   rc      - Planar radius of the base polygon on the tangent plane at the
%             north pole (same units as R).
%   l       - Number of polygon sides (integer, l >= 3).
%   phase   - In-plane rotation (radians) applied to the polygon.
%   R       - Sphere radius.
%
% Output:
%   P       - (l+1)-by-3 array of Cartesian points on the sphere; the last
%             point repeats the first to close the polygon.
%
% Notes:
%   1) The polygon is first built on the tangent plane at the north pole,
%      radially projected onto the sphere, and then rotated so that the
%      north pole aligns with center/||center|| using Rodrigues' formula.
%   2) For small rc/R, the polygon approximates a geodesic circle of angular
%      radius ≈ rc/R around 'center'.
%==========================================================================

    % Base polygon at the north pole (tangent plane z = R)
    ang = (0:l-1) * (2*pi/l) + phase;                     % Vertex angles
    P0  = [rc*cos(ang).', rc*sin(ang).', 0*ang.'] ...     % Planar l-gon
        + [0, 0, R];                                      % Shift to z = R

    % Radial projection of base polygon onto the sphere of radius R
    P0 = R * P0 ./ vecnorm(P0, 2, 2);

    % Rotation: align the north pole [0 0 1] with the target direction n
    n = center(:) / norm(center);                         % Unit target
    k = cross([0; 0; 1], n);                              % Rotation axis
    s = norm(k);                                          % |k| = sin(theta)
    c = dot([0; 0; 1], n);                                % cos(theta)

    if s < 1e-12
        Rm = eye(3);                                      % No rotation needed
    else
        % Rodrigues' rotation matrix from k (skew-symmetric K)
        K  = [   0,   -k(3),  k(2);
               k(3),    0,   -k(1);
              -k(2),  k(1),    0   ];
        Rm = eye(3) + K + K*K * ((1 - c) / (s^2));
    end

    % Apply rotation and close the polygon
    P = (Rm * P0.').';
    P(end + 1, :) = P(1, :);                              % Close polygon
end
