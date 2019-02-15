function [Fn, Ft] = decomp(F, D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Fn = sum(F .* D, 2) .* D;
Ft = F - Fn;
end

