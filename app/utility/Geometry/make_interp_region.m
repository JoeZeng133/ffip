function [x, y, z, dim, dx] = make_interp_region(p1, p2, varargin)
%MAKE_INTERP_REGION Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
addParameter(p, 'dim', [0 0 0]);
addParameter(p, 'dx', [0 0 0]);
parse(p, varargin{:});

if p.Results.dim(1) ~= 0
    dim = p.Results.dim(:)';
    dx = (p2(:) - p1(:)) ./ max(1, dim(:) - 1);
	dx = dx';
else
    if p.Results.dx(1) ~= 0
        dx = p.Results.dx;
        dim = round((p2(:) - p1(:)) ./ dx(:));
        dx = (p2(:) - p1(:)) ./ dim(:);
        dim = dim' + 1;
		dx = dx';
    else
        error('Either dx or dim should be provided')
    end
end

[x, y, z] = ndgrid(...
    linspace(p1(1), p2(1), dim(1)),...
    linspace(p1(2), p2(2), dim(2)),...
    linspace(p1(3), p2(3), dim(3)));
end

