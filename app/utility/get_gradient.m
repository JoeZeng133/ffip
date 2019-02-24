function [Fx, Fy, Fz] = get_gradient(f, varargin)
%GET_GRADIENT Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
addParameter(p, 'padding', 0);
addParameter(p, 'dx', [1 1 1]);
parse(p, varargin{:});

padded_f = padarray(f, [1 1 1], p.Results.padding, 'both');
dx = p.Results.dx;
[Fx, Fy, Fz] = gradient(padded_f, dx(1), dx(2), dx(3));
Fx = Fx(2:end - 1, 2:end - 1, 2:end - 1);
Fy = Fy(2:end - 1, 2:end - 1, 2:end - 1);
Fz = Fz(2:end - 1, 2:end - 1, 2:end - 1);

end

