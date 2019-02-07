function [poly, x, y] = make_disk(pos, radius, varargin)
%return a polygon of disk geometry centered at pos with radius
%   may specify the number of points on the perimiter
p = inputParser;
addParameter(p,'num_elements', 20);
parse(p, varargin{:});

t = linspace(0, 2 * pi, p.Results.num_elements);
x = pos(1) + cos(t(2:end)) * radius;
y = pos(2) + sin(t(2:end)) * radius;
poly = polyshape(x, y);


end

