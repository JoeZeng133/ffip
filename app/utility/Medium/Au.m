function [config, er, ur] = Au(varargin)
%Gold medium, intput unit in radian, valid for 200nm - 1000nm
%   Detailed explanation goes here
config = struct('er', 1.1431, 'sigma_e', 0, 'ur', 1, 'sigma_u', 0, 'poles', ...
    {{struct('type', 'Drude', 'freq', 1.3202e16/(2*pi), 'inv_relaxation', 1.0805e14/(2*pi)),...
    struct('type', 'CP', 'A', 0.26698, 'phi', -1.2371, 'Omega', 3.8711e15/(2*pi), 'Gamma', 4.4642e14/(2*pi)),...
    struct('type', 'CP', 'A', 3.0834, 'phi', -1.0968, 'Omega', 4.1684e15/(2*pi), 'Gamma', 2.3555e15/(2*pi))}});

if numel(varargin) == 0
    er = [];
    ur = [];
else
    w = varargin{1};
    er = config.er + ...
        drude(w, config.poles{1}.freq, config.poles{1}.inv_relaxation) + ...
        cp(w, config.poles{2}.A, config.poles{2}.phi, config.poles{2}.Omega, config.poles{2}.Gamma) + ...
        cp(w, config.poles{3}.A, config.poles{3}.phi, config.poles{3}.Omega, config.poles{3}.Gamma);
    ur = ones(size(w));
end

end

