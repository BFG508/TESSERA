function [tSec, qWm2] = loadHeatFlux(fileName, scaleFactor)
%==========================================================================
% loadHeatFlux: Reads thermal heat-flux from a spreadsheet-like file
%               [t(s), q(W/m^2)] and applies a scale factor. If the file
%               does not exist, a smooth synthetic profile is generated.
%
% Inputs:
%   fileName    - Path to a file readable by READMATRIX (e.g., .xlsx, .csv).
%   scaleFactor - Multiplicative scaling applied to q (default = 1).
%
% Outputs:
%   tSec - Time vector in seconds (column).
%   qWm2 - Heat-flux vector in W/m^2 (column), same length as tSec.
%
% Definition:
%   If the file exists, the function reads the first two numeric columns,
%   scales q by 'scaleFactor', sorts samples by time, and guarantees a
%   starting sample at t = 0 by prepending (0, q(1)) when necessary.
%   If the file is missing, it returns a synthetic profile (sum of two
%   Gaussian pulses) on [0, 250] s, also scaled by 'scaleFactor'.
%==========================================================================

    if nargin < 2 || isempty(scaleFactor)
        scaleFactor = 1;
    end

    if exist(fileName, 'file') == 2
        data = readmatrix(fileName);

        % Basic shape check (expect at least two columns)
        if isempty(data) || size(data,2) < 2
            warning('loadHeatFlux:InvalidData', ...
                'File "%s" has insufficient numeric data. Generating synthetic profile.', fileName);
            tSec = linspace(0,250,600).';
            qWm2 = scaleFactor * ( 5e4 * exp(-((tSec-60)./20).^2) ...
                                 + 2e4 * exp(-((tSec-140)./35).^2) );
            return;
        end

        % Extract and scale
        tSec = data(:,1);
        qWm2 = scaleFactor * data(:,2);

        % Cleaning: enforce columns, sort by time, and ensure t starts at 0
        [tSec, idx] = sort(tSec(:), 'ascend');
        qWm2 = qWm2(:);
        qWm2 = qWm2(idx);

        if tSec(1) > 0
            tSec = [0; tSec];
            qWm2 = [qWm2(1); qWm2];
        end

    else
        warning('loadHeatFlux:FileNotFound', ...
            'File "%s" not found. Generating synthetic heat flux.', fileName);

        tSec = linspace(0,250,600).';
        qWm2 = scaleFactor * ( 5e4 * exp(-((tSec-60)./20).^2) ...
                             + 2e4 * exp(-((tSec-140)./35).^2) );
    end
end
