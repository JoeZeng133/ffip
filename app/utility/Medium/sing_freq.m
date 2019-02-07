function [config, er, ur] = sing_freq(er, f, varargin)
% Match permittivity at one frequency, real(er) cannot be less than 0

e0 = 8.85418782e-12;
config = struct('er', real(er), 'sigma_e', -imag(er) * (2 * pi * f * e0), 'ur', 1, 'sigma_u', 0, 'poles', {});

if numel(varargin) == 0
    er = [];
    ur = [];
else
    w = varargin{1};
    er = real(er) - 1j * config.sigma_e ./ (w * e0);
    ur = ones(size(w));
end

end