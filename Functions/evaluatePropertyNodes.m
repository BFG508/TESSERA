function prop = evaluatePropertyNodes(layers, x, idx_s, idx_e, T, which)
%==========================================================================
% v: Evaluates a thermal property at mesh nodes by layer.
%
% Purpose:
%   Returns the nodal vector of a requested property ('k','cp','rho'),
%   evaluated per layer and as a function of temperature using simple models.
%
% Inputs:
%   layers : struct array with fields per layer:
%            - k_model/cp_model/rho_model : 'const' | 'lin' | 'poly2' (alias: 'poly')
%            - k_coeffs/cp_coeffs/rho_coeffs : [c0 c1 c2] (unused coeffs may be 0)
%   x      : (N×1) global node positions (m), concatenated across layers
%   idx_s  : (L×1) start node indices in x for each layer
%   idx_e  : (L×1) end   node indices in x for each layer
%   T      : (N×1) nodal temperatures (K) at current time level
%   which  : char, one of 'k','cp','rho' (property to evaluate)
%
% Output:
%   prop   : (N×1) nodal property vector (SI units), clamped to ≥ eps
%
% Models:
%   const : p(T) = c0
%   lin   : p(T) = c0 + c1*T
%   poly2 : p(T) = c0 + c1*T + c2*T^2   (also accepts alias 'poly')
%==========================================================================

    % Basic checks (lightweight)
    if numel(T) ~= numel(x)
        error('eval_property_nodes:SizeMismatch', ...
              'Length of T (%d) must match length of x (%d).', numel(T), numel(x));
    end

    prop = zeros(numel(x),1);

    for i = 1:numel(layers)
        ids = idx_s(i):idx_e(i);
        Ti  = T(ids);

        % Select property model and coeffs for this layer
        switch lower(which)
            case 'k'
                model = lower(layers(i).k_model);
                c     = layers(i).k_coeffs(:).';
            case 'cp'
                model = lower(layers(i).cp_model);
                c     = layers(i).cp_coeffs(:).';
            case 'rho'
                model = lower(layers(i).rho_model);
                c     = layers(i).rho_coeffs(:).';
            otherwise
                error('eval_property_nodes:UnknownProperty', ...
                      'Unknown property "%s". Expected ''k'', ''cp'', or ''rho''.', which);
        end

        % Evaluate model (vectorized over nodes in this layer)
        switch model
            case 'const'
                % p(T) = c0
                prop(ids) = c(1);

            case 'lin'
                % p(T) = c0 + c1*T
                prop(ids) = c(1) + c(2)*Ti;

            case {'poly2','poly'}
                % p(T) = c0 + c1*T + c2*T^2
                prop(ids) = c(1) + c(2)*Ti + c(3)*Ti.^2;

            otherwise
                error('eval_property_nodes:UnsupportedModel', ...
                      'Model "%s" not supported for property "%s".', model, which);
        end

        % Prevent non-physical / zero values
        prop(ids) = max(prop(ids), eps);
    end
end
