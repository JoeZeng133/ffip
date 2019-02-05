function [config, er, ur] = Air(varargin)
%air medium, intput unit in radian
%   Detailed explanation goes here
config = struct('er', 1, 'sigma_e', 0, 'ur', 1, 'sigma_u', 0, 'poles', {{}});

if numel(varargin) == 0
    er = [];
    ur = [];
else
    er = ones(size(varargin{1}));
    ur = ones(size(varargin{1}));
end

end

