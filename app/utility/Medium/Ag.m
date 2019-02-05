function [config, er, ur] = Ag(varargin)
% Silver, valid for 200nm - 1000nm

config = struct('er', 1.4447, 'sigma_e', 0, 'ur', 1, 'sigma_u', 0, 'poles', ...
    {{struct('type', 'Drude', 'freq', 1.3280e16/(2*pi), 'inv_relaxation', 9.1269e13/(2*pi)),...
    struct('type', 'CP', 'A', -1.5951, 'phi', 3.1288, 'Omega', 8.2749e15/(2*pi), 'Gamma', 5.177e15/(2*pi)),...
    struct('type', 'CP', 'A', 0.25261, 'phi', -1.5066, 'Omega', 6.1998e15/(2*pi), 'Gamma', 5.4126e14/(2*pi))}});

if numel(varargin) == 0
    er = [];
    ur = [];
else
    w = varargin{1};
    er = config.er+ ...
        drude(w, config.poles{1}.freq, config.poles{1}.inv_relaxation) + ...
        cp(w, config.poles{2}.A, config.poles{2}.phi, config.poles{2}.Omega, config.poles{2}.Gamma) + ...
        cp(w, config.poles{3}.A, config.poles{3}.phi, config.poles{3}.Omega, config.poles{3}.Gamma);
    ur = ones(size(w));
end

end