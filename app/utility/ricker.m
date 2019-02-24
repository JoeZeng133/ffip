function [res] = ricker(t, fp, d)
%RICKER Summary of this function goes here
%   Detailed explanation goes here
res = (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
end

