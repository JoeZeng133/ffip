function [x, y, z, dim] = make_interp_region(p1, p2, dx)
%MAKE_INTERP_REGION Summary of this function goes here
%   Detailed explanation goes here
dim = round((p2 - p1) / dx);
[x, y, z] = ndgrid(...
    linspace(p1(1), p2(1), dim(1) + 1),...
    linspace(p1(2), p2(2), dim(2) + 1),...
    linspace(p1(3), p2(3), dim(3) + 1));
dim = dim + 1;
end

