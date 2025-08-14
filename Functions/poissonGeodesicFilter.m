function [keep_pts, keep_lados, keep_theta] = poissonGeodesicFilter(pts, lados_vec, theta_vec, safety_fac)
%==========================================================================
% poissonGeodesicFilter: Poisson-disk filtering on the sphere using a
%                        lat/lon grid for near-neighbor acceleration.
%
% Inputs:
%   pts         - N-by-3 Cartesian candidate points on (or near) a sphere.
%   lados_vec   - N-by-1 auxiliary data vector associated with each point
%                 (e.g., number of sides); filtered to accepted points.
%   theta_vec   - N-by-1 angular radius per point (radians).
%   safety_fac  - Safety factor for exclusion distance (>= 1 recommended).
%
% Outputs:
%   keep_pts    - M-by-3 accepted points (M <= N).
%   keep_lados  - M-by-1 filtered entries from lados_vec.
%   keep_theta  - M-by-1 filtered entries from theta_vec.
%
% Method:
%   Greedy acceptance: points are processed in input order. Each accepted
%   point excludes future points within an angular distance less than
%   (theta_i + theta_j) * safety_fac. To avoid O(N^2) checks, candidates are
%   binned into a latitude/longitude grid whose cell size is tied to the
%   minimum angular radius. Only neighbors in adjacent cells are compared.
%
% Notes:
%   - Points are normalized to unit vectors before angular computations.
%   - Azimuth indices wrap around; latitude indices are clamped.
%   - The dot product is clamped to [-1, 1] before acos for robustness.
%   - Expected complexity is ~O(N), assuming roughly uniform distribution.
%==========================================================================

    % Basic sizes and unit vectors
    N  = size(pts, 1);
    R  = vecnorm(pts, 2, 2);                 % Radius per point
    U  = pts ./ R;                           % Normalize to unit vectors

    % Spherical coordinates
    phi = acos(U(:, 3));                     % Polar angle in [0, pi]
    lam = atan2(U(:, 2), U(:, 1));           % Azimuth in [-pi, pi)

    % Angular cell size: use minimum exclusion radius (scaled) for coverage
    cell_ang = max(1e-6, min(theta_vec) * safety_fac);

    % Grid dimensions
    nPhi = max(1,  floor(max(phi) / cell_ang) + 1);       % # lat cells
    nLam = max(12, floor(2*pi     / cell_ang));           % >= 12 azimuth

    % Cell indices per point (clamped in latitude, wrapped in azimuth)
    dPhi = max(phi) / max(nPhi - 1, 1) + eps;             % Avoid division by zero
    dLam = 2*pi / nLam;

    iPhi = min(nPhi - 1, max(0, floor(phi / dPhi)));      % [0, nPhi-1]
    iLam = floor((lam + pi) / dLam);                       % [0, nLam]
    iLam(iLam == nLam) = nLam - 1;                         % Wrap to [0, nLam-1]

    % Hash buckets: map (iPhi, iLam) -> list of accepted indices
    buckets = cell(nPhi, nLam);

    % Greedy acceptance state
    keep            = false(N, 1);
    accepted_U      = zeros(0, 3);                         % Unit vectors of accepted points
    accepted_theta  = zeros(0, 1);                         % Angular radii of accepted points

    % Main loop over candidates
    for i = 1:N
        ip = iPhi(i) + 1;                                  % 1-based latitude cell
        il = iLam(i) + 1;                                  % 1-based azimuth cell
        theta_i = theta_vec(i);
        u = U(i, :);

        % Neighborhood (in cells) based on current angular radius
        rPhi = ceil(theta_i / dPhi);
        rLam = ceil(theta_i / dLam);

        ok = true;
        for dp = -rPhi:rPhi
            ip2 = ip + dp;
            if ip2 < 1 || ip2 > nPhi, continue; end

            for dl = -rLam:rLam
                il2 = il + dl;
                % Azimuth wrap-around
                if il2 < 1
                    il2 = il2 + nLam;
                elseif il2 > nLam
                    il2 = il2 - nLam;
                end

                idxList = buckets{ip2, il2};
                if isempty(idxList), continue; end

                % Exact angular checks within the neighboring bucket
                for jj = 1:numel(idxList)
                    j = idxList(jj);
                    % Clamp dot product for numerical robustness
                    ang = acos(max(-1, min(1, dot(u, accepted_U(j, :)))));
                    if ang < (theta_i + accepted_theta(j)) * safety_fac
                        ok = false; break;
                    end
                end
                if ~ok, break; end
            end
            if ~ok, break; end
        end

        % Accept current point if all checks passed
        if ok
            keep(i) = true;

            % Register acceptance
            accepted_U(end+1, :)     = u;
            accepted_theta(end+1, 1) = theta_i;
            buckets{ip, il} = [buckets{ip, il}, size(accepted_U, 1)];
        end
    end

    % Filter outputs to accepted points
    keep_pts   = pts(keep, :);
    keep_lados = lados_vec(keep);
    keep_theta = theta_vec(keep);
end