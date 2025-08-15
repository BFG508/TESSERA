function TPS = betaTPS()
%==========================================================================
% betaTPS: Build a ready-to-run config for transientConduction_multilayer
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

    % --- Layer 1: ceramic outer tile (low-k) ------------------------------
    layer1 = struct( ...
        'name','Tile', ...
        'L',   0.006, ...                % thickness [m]
        'Nx',  19,   ...                 % nodes
        'k_model',  'lin',  'k_coeffs',  [0.06, 1.0e-4, 0], ...  % W/m-K
        'cp_model', 'lin',  'cp_coeffs', [700,  0.20,   0], ...  % J/kg-K
        'rho_model','const','rho_coeffs',[180,  0,      0] ...   % kg/m^3
    );

    % --- Layer 2: fibrous insulation -------------------------------------
    layer2 = struct( ...
        'name','Insulation', ...
        'L',   0.012, ...
        'Nx',  25, ...
        'k_model',  'const','k_coeffs',  [0.040, 0, 0], ...
        'cp_model', 'lin',  'cp_coeffs', [900,  0.15, 0], ...
        'rho_model','const','rho_coeffs',[100,  0,    0] ...
    );

    % --- Layer 3: adhesive ------------------------------------------------
    layer3 = struct( ...
        'name','Adhesive', ...
        'L',   0.0015, ...
        'Nx',  5, ...
        'k_model',  'const','k_coeffs',  [0.35, 0, 0], ...
        'cp_model', 'const','cp_coeffs', [1000, 0, 0], ...
        'rho_model','const','rho_coeffs',[1200, 0, 0] ...
    );

    % --- Layer 4: aluminum backshell -------------------------------------
    layer4 = struct( ...
        'name','Aluminum', ...
        'L',   0.004, ...
        'Nx',  9, ...
        'k_model',  'lin',  'k_coeffs',  [120.0, -1.0e-2, 0], ... % mild negative slope
        'cp_model', 'const','cp_coeffs', [900,   0,      0], ...
        'rho_model','const','rho_coeffs',[2700,  0,      0] ...
    );

    % --- Time vector ------------------------------------------------------
    % 0 to 400 s with ~0.2 s resolution
    t_vec = linspace(0, 400, 2001)';  % [s]

    % --- Assemble TPS -----------------------------------------------------
    TPS = struct();
    TPS.layers = [layer1, layer2, layer3, layer4];
    TPS.time   = struct('t', t_vec);

    % --- Hot-side BC: prescribed heat flux q_in(t) [W/m^2] ---------------
    % Two Gaussian pulses + baseline to mimic a multi-peak heating profile
    TPS.hotBC = struct( ...
        'type','flux', ...
        'qfun', @(tt) (7.0e4*exp(-((tt-150)/40).^2) + ...
                       3.0e4*exp(-((tt-280)/60).^2) + ...
                       5.0e3) ...
    );

    % --- Cold-side BC: insulated (Neumann 0) ------------------------------
    TPS.coldBC = struct('type','insulated');

    % --- Energy area (optional) ------------------------------------------
    TPS.area = pi*(0.25/2)^2;  % m^2, e.g., 0.25 m diameter disk (demo)

    % --- Store temperature history ---------------------------------------
    TPS.save_history = true;

    % --- Initial temperature ---------------------------------------------
    TPS.T_init = 290.0;  % K
end
