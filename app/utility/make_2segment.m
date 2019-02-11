function [res] = make_2segment(a, b, len)
%MAKE_2SEGMENT Summary of this function goes here
%   Detailed explanation goes here
if len > 2
    res = b * ones([len 1]);
    res(1) = a;
    res(end) = a;
else
    if (len == 2)
        res = [a; a];
    end
    
    if (len == 1)
        res = b;
    end
end
end

