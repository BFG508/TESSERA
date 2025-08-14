function pts = fibonacciPoints(N, R, phi_max)
%==========================================================================
% fibonacci_cap_points: Generates quasi-uniform Fibonacci sampling of points
%                       on a spherical cap of half-angle phi_max.
%
% Inputs:
%   N        - Number of points (nonnegative integer).
%   R        - Sphere radius.
%   phi_max  - Polar half-angle of the spherical cap in radians (0 <= phi_max <= pi).
%
% Output:
%   pts      - N-by-3 array of Cartesian points on the sphere of radius R
%              with 0 <= phi <= phi_max.
%
% Notes:
%   Uses a Fibonacci spiral with the golden angle increment. Points are sampled
%   uniformly in z over [cos(phi_max), 1] and uniformly in azimuth via the
%   golden-angle rotation.
%==========================================================================

    % Handle empty/degenerate request
    if N < 1
        pts = zeros(0, 3);  % Return an empty N×3 array
        return;
    end

    % Golden angle for Fibonacci spiral placement
    ga = goldenAngle;

    % Lower bound in z for the cap and preallocate output
    zmin = cos(phi_max);     % Cap limit at polar angle phi_max
    pts  = zeros(N, 3);      % Preallocate N×3 output

    % Generate N points on the spherical cap
    for i = 0:N-1
        % z uniformly distributed over [zmin, 1]
        z = zmin + (1 - zmin) * (i + 0.5) / N;  % Midpoint rule in z

        % Radius of the projection on the XY-plane
        r_xy = sqrt(max(0, 1 - z^2));           % Clamp for round-off

        % Azimuth using the golden-angle increment
        theta = ga * i;

        % Cartesian coordinates on the unit sphere, then scale by R
        x = r_xy * cos(theta);
        y = r_xy * sin(theta);
        pts(i + 1, :) = R * [x, y, z];
    end

end
