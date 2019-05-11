function [x, y, z, f, v] = read_dft(filename, groupname)
%READ_DFT Summary of this function goes here
%   Detailed explanation goes here
x = h5read(filename, [groupname, '/x']);
y = h5read(filename, [groupname, '/y']);
z = h5read(filename, [groupname, '/z']);
f = h5read(filename, [groupname, '/f']);
v = h5read(filename, [groupname, '/real']) ...
    + 1j * h5read(filename, [groupname, '/imag']);

end

