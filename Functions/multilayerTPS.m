function out = multilayerTPS(cfg)
%==========================================================================
% multilayerTPS: Transient 1D multilayer heat solver (TPS)
%
% Purpose:
%   Solve transient 1D heat conduction through a stack of TPS layers using
%   a Crank–Nicolson control-volume formulation with temperature-dependent
%   properties per layer.
%
% Inputs:
%   cfg                : struct with simulation setup:
%     .layers(i)       : struct per layer with fields
%                        - name        : char
%                        - L           : thickness [m]
%                        - Nx          : number of nodes (>=2)
%                        - k_model     : 'const' | 'lin' | 'poly2' (alias: 'poly')
%                        - k_coeffs    : [c0 c1 c2] (unused coeffs may be 0)
%                        - cp_model    : idem
%                        - cp_coeffs   : idem
%                        - rho_model   : idem
%                        - rho_coeffs  : idem
%     .time.t          : (Nt×1) time vector [s], strictly increasing
%     .hotBC           : struct with hot-side boundary condition
%                        type='flux'  → .qfun(t) [W/m^2]
%                        type='robin' → .h [W/m^2K], .Tinf [K]
%     .coldBC          : struct with cold-side BC
%                        'robin' (.h,.Tinf) | 'insulated' | 'dirichlet' (.Tfix)
%     .area            : surface area [m^2] for energy bookkeeping (optional)
%     .save_history    : true/false (default true)
%     .T_init          : scalar or (Nnodes×1) initial temperature [K]
%
% Outputs:
%   out                : struct with results
%     .x               : (N×1) global node positions [m]
%     .t               : (Nt×1) time vector [s]
%     .T               : (Nt×N) temperature history [K] (if save_history==true)
%     .T_end           : (N×1) final temperature profile [K]
%     .T_hot           : (Nt×1) hot-side temperature [K]
%     .T_cold          : (Nt×1) cold-side temperature [K]
%     .T_ifaces        : (Nt×(nlayers-1)) interface temperatures [K]
%     .meta            : layer info and index maps
%     .energy          : struct with Ein, Eout_cold, U, residual (if applicable)
%
% Method:
%   Crank–Nicolson (second order in time). Control volumes per node with
%   thermal capacities C_i = (rho*cp)_i * Δx_i and face conductances
%   via harmonic average of k and geometric distances. Material properties
%   are evaluated at T^n (semi-implicit linearization).
%
% Notes:
%   - Nonlinear radiation is not included in this base version (can be added
%     via εσ linearization). Use hotBC.type='flux' for net surface flux.
%   - Requires the helper: eval_property_nodes(...)
%==========================================================================

    %% --------------------- Basic validation ---------------------
    assert(isfield(cfg,'layers') && ~isempty(cfg.layers), ...
        'cfg.layers is required and must be non-empty.');
    assert(isfield(cfg,'time') && isfield(cfg.time,'t'), ...
        'cfg.time.t is required.');
    t  = cfg.time.t(:);
    Nt = numel(t);
    assert(Nt >= 2, 'cfg.time.t must contain at least two time points and be increasing.');

    % Hot-side BC
    hot = cfg.hotBC;
    cold = cfg.coldBC;
    if strcmpi(hot.type,'flux')
        assert(isfield(hot,'qfun') && isa(hot.qfun,'function_handle'), ...
              'For hotBC.type=''flux'', provide hotBC.qfun(t) as a function handle.');
    elseif strcmpi(hot.type,'robin')
        assert(all(isfield(hot,{'h','Tinf'})), ...
              'For hotBC.type=''robin'', provide hotBC.h and hotBC.Tinf.');
    else
        error('Invalid hotBC.type. Use ''flux'' or ''robin''.');
    end

    % Cold-side BC
    switch lower(cold.type)
        case 'robin'
            assert(all(isfield(cold,{'h','Tinf'})), ...
                  'For coldBC.type=''robin'', provide coldBC.h and coldBC.Tinf.');
        case 'insulated'
            % OK (Neumann 0)
        case 'dirichlet'
            assert(isfield(cold,'Tfix'), ...
                  'For coldBC.type=''dirichlet'', provide coldBC.Tfix.');
        otherwise
            error('Invalid coldBC.type. Use ''robin'', ''insulated'', or ''dirichlet''.');
    end

    Aref = 1.0;
    if isfield(cfg,'area') && ~isempty(cfg.area)
        Aref = cfg.area;
    end

    save_hist = true;
    if isfield(cfg,'save_history'), save_hist = logical(cfg.save_history); end

    %% --------------------- Build multilayer mesh ---------------------
    Ls  = arrayfun(@(L) L.L,  cfg.layers);
    Nxs = arrayfun(@(L) L.Nx, cfg.layers);
    assert(all(Ls > 0) && all(Nxs >= 2), ...
          'Each layer must have L>0 and Nx>=2.');

    % Local nodes per layer in [0, L_layer]
    x_local = cell(numel(cfg.layers),1);
    for i = 1:numel(cfg.layers)
        x_local{i} = linspace(0, cfg.layers(i).L, cfg.layers(i).Nx);
    end

    % Shift and concatenate, removing duplicated interface nodes
    offsets = [0, cumsum(Ls(1:end-1))];
    for i = 1:numel(cfg.layers)
        x_local{i} = x_local{i} + offsets(i);
    end
    x = x_local{1};
    idx_layer_start = zeros(numel(cfg.layers),1);
    idx_layer_end   = zeros(numel(cfg.layers),1);
    idx_layer_start(1) = 1;
    for i = 2:numel(cfg.layers)
        x = [x, x_local{i}(2:end)]; %#ok<AGROW>
        idx_layer_start(i) = numel(x_local{i-1}) - 1 + idx_layer_start(i-1);
    end
    idx_layer_end = [idx_layer_start(2:end)-1; numel(x)];
    x  = x(:);
    Nn = numel(x);

    % Control-volume lengths Δx_i (center-to-midface distances)
    dx = zeros(Nn,1);
    dx(1)   = (x(2) - x(1))/2;
    dx(end) = (x(end) - x(end-1))/2;
    for i = 2:Nn-1
        dx(i) = (x(i+1) - x(i-1))/2;
    end

    % Interface node indices (last node of each layer except the last layer)
    iface_idx = idx_layer_end(1:end-1);

    %% --------------------- Initial temperature ---------------------
    if isfield(cfg,'T_init') && ~isempty(cfg.T_init)
        if isscalar(cfg.T_init)
            T = cfg.T_init * ones(Nn,1);
        else
            assert(numel(cfg.T_init) == Nn, 'cfg.T_init length must match number of nodes.');
            T = cfg.T_init(:);
        end
    else
        % Default: average of available ambient values
        T0 = 300;
        if strcmpi(hot.type,'robin') && isfield(hot,'Tinf'),   T0 = hot.Tinf;         end
        if strcmpi(cold.type,'robin') && isfield(cold,'Tinf'), T0 = 0.5*(T0+cold.Tinf); end
        if strcmpi(cold.type,'dirichlet'),                     T0 = 0.5*(T0+cold.Tfix); end
        T = T0 * ones(Nn,1);
    end

    %% --------------------- History buffers ---------------------
    if save_hist
        T_hist = zeros(Nt, Nn); T_hist(1,:) = T.';
    else
        T_hist = [];
    end
    T_hot  = zeros(Nt,1); T_hot(1)  = T(1);
    T_cold = zeros(Nt,1); T_cold(1) = T(end);
    T_if   = zeros(Nt, numel(iface_idx));
    if ~isempty(iface_idx), T_if(1,:) = T(iface_idx); end

    % Energy bookkeeping (only meaningful for hotBC='flux' and known area)
    do_energy   = (strcmpi(hot.type,'flux') && Aref > 0);
    E_in        = 0;
    E_out_cold  = 0;

    %% ===================== Time-stepping (Crank–Nicolson) =====================
    for n = 1:Nt-1
        dt = t(n+1) - t(n);
        if dt <= 0
            error('Time vector must be strictly increasing.');
        end

        % Properties at T^n (semi-implicit)
        k_node   = evaluatePropertyNodes(cfg.layers, x, idx_layer_start, idx_layer_end, T, 'k');
        cp_node  = evaluatePropertyNodes(cfg.layers, x, idx_layer_start, idx_layer_end, T, 'cp');
        rho_node = evaluatePropertyNodes(cfg.layers, x, idx_layer_start, idx_layer_end, T, 'rho');

        % Nodal capacities per unit area
        C = rho_node .* cp_node .* dx;

        % Face conductances per unit area using harmonic average of k
        G = zeros(Nn-1,1);
        for i = 1:Nn-1
            kL = k_node(i);
            kR = k_node(i+1);
            dL = (x(i+1)-x(i))/2;  % node-to-face distance on the left
            dR = dL;               % same on the right for uniform subsegments
            Rth = dL/max(kL,eps) + dR/max(kR,eps);
            G(i) = 1.0 / Rth;
        end

        % Assemble tridiagonal matrices A*T^{n+1} = B*T^n + rhs  (CN scheme)
        A   = spalloc(Nn, Nn, 3*Nn);
        B   = spalloc(Nn, Nn, 3*Nn);
        rhs = zeros(Nn,1);

        % Interior nodes
        for i = 2:Nn-1
            GiL = G(i-1); GiR = G(i);
            A(i,i-1) = -0.5*dt*GiL;
            A(i,i)   =  C(i) + 0.5*dt*(GiL + GiR);
            A(i,i+1) = -0.5*dt*GiR;

            B(i,i-1) =  0.5*dt*GiL;
            B(i,i)   =  C(i) - 0.5*dt*(GiL + GiR);
            B(i,i+1) =  0.5*dt*GiR;
        end

        % ---- Hot boundary (i = 1) ----
        GiR = G(1);
        A(1,1) = C(1) + 0.5*dt*GiR;
        A(1,2) =       -0.5*dt*GiR;
        B(1,1) = C(1) - 0.5*dt*GiR;
        B(1,2) =        0.5*dt*GiR;

        if strcmpi(hot.type,'flux')
            qn  = hot.qfun(t(n));
            qn1 = hot.qfun(t(n+1));
            rhs(1) = rhs(1) + 0.5*dt*(qn + qn1);  % trapezoidal in time
        elseif strcmpi(hot.type,'robin')
            h = hot.h; Tinf = hot.Tinf;
            % Implicit in T, explicit constant part in rhs
            A(1,1) = A(1,1) + 0.5*dt*h;
            B(1,1) = B(1,1) - 0.5*dt*h;
            rhs(1) = rhs(1) + 0.5*dt*h*(Tinf + Tinf);
        end

        % ---- Cold boundary (i = Nn) ----
        GiL = G(end);
        switch lower(cold.type)
            case 'robin'
                h = cold.h; Tinf = cold.Tinf;
                A(Nn,Nn-1) = -0.5*dt*GiL;
                A(Nn,Nn)   =  C(Nn) + 0.5*dt*(GiL + h);
                B(Nn,Nn-1) =  0.5*dt*GiL;
                B(Nn,Nn)   =  C(Nn) - 0.5*dt*(GiL + h);
                rhs(Nn)    =  rhs(Nn) + 0.5*dt*h*(Tinf + Tinf);
            case 'insulated' % Neumann 0
                A(Nn,Nn-1) = -0.5*dt*GiL;
                A(Nn,Nn)   =  C(Nn) + 0.5*dt*GiL;
                B(Nn,Nn-1) =  0.5*dt*GiL;
                B(Nn,Nn)   =  C(Nn) - 0.5*dt*GiL;
            case 'dirichlet'
                % Enforce T_N^{n+1} = Tfix
                A(Nn,:) = 0; A(Nn,Nn) = 1;
                B(Nn,:) = 0;
                rhs(Nn) = cold.Tfix;
        end

        % Solve linear system
        b    = B*T + rhs;
        Tnew = A \ b;

        % Energy bookkeeping (optional)
        if do_energy
            qn   = hot.qfun(t(n));
            qn1  = hot.qfun(t(n+1));
            E_in_step = 0.5*(qn + qn1) * Aref * dt;   % energy IN at hot side
            E_in = E_in + E_in_step;

            if strcmpi(cold.type,'robin')
                % Heat leaving at cold side (>0 if Tsurf > Tinf)
                Tsurf_cold = Tnew(end);
                q_out = cold.h * (Tsurf_cold - cold.Tinf);
                E_out_cold = E_out_cold + max(0, q_out) * Aref * dt;
            end
        end

        % Advance and store history
        T = Tnew;
        if save_hist, T_hist(n+1,:) = T.'; end
        T_hot(n+1)  = T(1);
        T_cold(n+1) = T(end);
        if ~isempty(iface_idx), T_if(n+1,:) = T(iface_idx); end
    end

    %% --------------------------- Output ---------------------------
    out = struct();
    out.x        = x;
    out.t        = t;
    out.T        = T_hist;
    out.T_end    = T;
    out.T_hot    = T_hot;
    out.T_cold   = T_cold;
    out.T_ifaces = T_if;
    out.meta.layers          = cfg.layers;
    out.meta.idx_layer_start = idx_layer_start;
    out.meta.idx_layer_end   = idx_layer_end;

    if do_energy
        % Internal energy at final time (reference at 0 K without loss of generality)
        % Uses last-step rho_node, cp_node in scope
        U = sum((rho_node .* cp_node .* dx) .* T);
        out.energy = struct('Ein', E_in, 'Eout_cold', E_out_cold, ...
                            'U', U, 'residual', E_in - E_out_cold - U);
    else
        out.energy = struct();
    end
end