function TPS = alfaTPS()
%==========================================================================
% alfaTPS: Build a ready-to-run config for transientConduction_multilayer
%
% Purpose:
%   Provide a minimal, self-contained configuration (layers, time vector,
%   and boundary conditions) to test transientConduction_multilayer.m.
%
% Output:
%   TPS : struct with fields
%         .layers(i)  : layer definition (see below)
%         .time.t     : (Nt×1) time vector [s], strictly increasing
%         .hotBC      : hot-side boundary condition (flux or Robin)
%         .coldBC     : cold-side boundary condition (Robin/insulated/Dirichlet)
%         .area       : surface area [m^2] used for energy bookkeeping
%         .save_history : logical flag to store T(x,t)
%         .T_init     : initial temperature (scalar or nodal vector) [K]
%
% Layer fields (per layer):
%   name, L [m], Nx (>=2), k_model, k_coeffs [c0 c1 c2],
%   cp_model, cp_coeffs [c0 c1 c2], rho_model, rho_coeffs [c0 c1 c2].
%   Models: 'const' | 'lin' | 'poly2' (alias 'poly' accepted by the solver).
%==========================================================================

    % --- Layer 1: tile (insulation) --------------------------------------
    layer1 = struct( ...
        'name','Tile', ...
        'L',   0.010, ...                 % thickness [m]
        'Nx',  21,    ...                 % nodes (>=2)
        'k_model',  'poly2', 'k_coeffs',  [0.05, 1e-4, 0], ...   % W/m-K
        'cp_model', 'poly2', 'cp_coeffs', [900,  0.1,  0], ...   % J/kg-K
        'rho_model','const', 'rho_coeffs',[250,  0,    0] ...    % kg/m^3
    );

    % --- Layer 2: adhesive ------------------------------------------------
    layer2 = struct( ...
        'name','Adhesive', ...
        'L',   0.002, ...
        'Nx',  7, ...
        'k_model',  'const', 'k_coeffs',  [0.30, 0, 0], ...
        'cp_model', 'const', 'cp_coeffs', [1000, 0, 0], ...
        'rho_model','const', 'rho_coeffs',[1100, 0, 0] ...
    );

    % --- Layer 3: backshell ----------------------------------------------
    layer3 = struct( ...
        'name','Backshell', ...
        'L',   0.015, ...
        'Nx',  21, ...
        'k_model',  'poly2', 'k_coeffs',  [16.0, 0.01, 0], ...
        'cp_model', 'const', 'cp_coeffs', [880,  0,    0], ...
        'rho_model','const', 'rho_coeffs',[2700, 0,    0] ...
    );

    % --- Time vector ------------------------------------------------------
    % 0 to 150 s with fine resolution (10,000 points)
    t_vec = linspace(0, 150, 1e4)';   % [s], strictly increasing

    % --- Assemble TPS -----------------------------------------------------
    TPS = struct();
    TPS.layers       = [layer1, layer2, layer3];
    TPS.time         = struct('t', t_vec);

    % --- Hot-side BC: prescribed heat flux q_in(t) [W/m^2] ---------------
    TPS.hotBC = struct( ...
        'type','flux', ...
        'qfun', @(tt) 5e4 * exp(-((tt-60)/20).^2) ...   % demo Gaussian pulse
    );

    % --- Energy area (for bookkeeping, optional) -------------------------
    TPS.area = pi*(0.4/2)^2;  % m^2, e.g., 0.4 m diameter disk (demo)

    % --- Store temperature history ---------------------------------------
    TPS.save_history = true;

    % --- Initial temperature ---------------------------------------------
    TPS.T_init = 303.15;      % K (scalar; solver expands to all nodes)

    % --- Cold-side BC: convective Robin ----------------------------------
    TPS.coldBC = struct( ...
        'type',  'robin', ...
        'h',     10, ...        % W/m^2-K
        'Tinf',  TPS.T_init ...     % K
    );
end
