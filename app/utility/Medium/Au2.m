function [config, er, ur] = Au2(varargin)
% Au medium valid for 600nm-1000nm

config = struct('er', 5.9673, 'sigma_e', 0, 'ur', 1, 'sigma_u', 0, 'poles', ...
    {{struct('type', 'Drude', 'freq', 2113.6e12, 'inv_relaxation', 15.92e12),...
    struct('type', 'Lorentz', 'rel_perm', 1.09, 'freq', 650.07e12, 'damp', 104.86e12)}});

if numel(varargin) == 0
    er = [];
    ur = [];
else
    w = varargin{1};
    er = config.er + drude(w, config.poles{1}.freq, config.poles{1}.inv_relaxation) + ...
        lorentz(w, config.poles{2}.rel_perm, config.poles{2}.freq, config.poles{2}.damp);
    ur = ones(size(w));
end

end