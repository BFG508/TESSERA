function A = sphericalPolygonArea(P, R)
%==========================================================================
% sphericalPolygonArea: Computes the area of a (closed) spherical polygon
%                       on a sphere of radius R.
%
% Inputs:
%   P   - N-by-3 array of polygon vertices (Cartesian coordinates). If the
%         polygon is not closed (P(1,:) ~= P(end,:)), it will be closed
%         internally by repeating the first vertex at the end.
%   R   - Sphere radius.
%
% Output:
%   A   - Spherical area (same units as R^2).
%
% Method:
%   The polygon is triangulated into a fan around an approximate center c0.
%   Each spherical triangle (with unit vectors a, b, c) has area on the unit
%   sphere given by the robust formula:
%       A_tri_unit = 2*atan2( |a · (b × c)| , 1 + a·b + b·c + c·a )
%   The final area is A = R^2 * sum(A_tri_unit).
%
% Notes:
%   - Vertices are normalized to unit vectors before area evaluation.
%   - The fan center c0 is taken as the normalized mean of vertices; for
%     highly concave polygons or polygons spanning > π steradians, a more
%     careful choice of fan center may be required.
%==========================================================================

    % Ensure the polygon is closed
    if ~isequal(P(1, :), P(end, :))
        P = [P; P(1, :)];
    end

    % Normalize vertices to unit vectors on the sphere
    U = P(1:end-1, :);
    U = U ./ vecnorm(U, 2, 2);           % Unit directions per vertex

    % Choose a fan center and normalize it
    c0 = mean(U, 1);
    c0 = c0 / norm(c0);

    % Accumulate unit-sphere area via triangle fan
    A_unit = 0;
    for k = 1:size(U, 1) - 1
        a = U(k,   :);
        b = U(k+1, :);
        c = c0;

        % Robust spherical triangle area on the unit sphere
        numer = abs(dot(a, cross(b, c)));
        denom = 1 + dot(a, b) + dot(b, c) + dot(c, a);
        A_unit = A_unit + 2 * atan2(numer, denom);
    end

    % Scale by R^2 to obtain area on sphere of radius R
    A = R^2 * A_unit;
end
