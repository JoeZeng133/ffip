function [x, y, z] = read_geom(filename, groupname)
%READ_GEOM Summary of this function goes here
%   Detailed explanation goes here
x = h5read(filename, [groupname, '/x']);
y = h5read(filename, [groupname, '/y']);
z = h5read(filename, [groupname, '/z']);
end

