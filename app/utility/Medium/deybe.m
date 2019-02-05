function [res] = deybe(w, rel_e, tau)
%cp pole, tau in sec
%   Detailed explanation goes here
res = rel_e ./ (1 + 1j * w * tau);
end

