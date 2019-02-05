function [config, er, ur] = fic1(varargin)
% fictious medium
config = struct('er', 1, 'sigma_e', 0, 'ur', 1, 'sigma_u', 0, 'poles', ...
{{struct('type', 'Lorentz', 'rel_perm', 3, 'freq', 3e14, 'damp', 0.5e14),...
struct('type', 'Lorentz', 'rel_perm', 3, 'freq', 5e14, 'damp', 1e14)}});

if numel(varargin) == 0
    er = [];
    ur = [];
else
    w = varargin{1};
    er = config.er + ...
        lorentz(w, config.poles{1}.rel_perm, config.poles{1}.freq, config.poles{1}.damp) + ...
        lorentz(w, config.poles{2}.rel_perm, config.poles{2}.freq, config.poles{2}.damp);
    ur = ones(size(w));
end

end