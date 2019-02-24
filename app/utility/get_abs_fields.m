function [res] = get_abs_fields(F)
%GET_ABS_FIELDS Summary of this function goes here
%   Detailed explanation goes here
if size(F, 2) ~= 3
    error('Invalid Fields');
end

res = sqrt(sum(abs(F).^2, 2));
end

