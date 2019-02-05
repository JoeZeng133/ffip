function [res] = lorentz(w, rel_e, fp, delta)
%lorentz pole, w in radians, fp delta in Hz
%   Detailed explanation goes here
wp = 2 * pi * fp;
delta = 2 * pi * delta;
res = rel_e ./ (1 + 1j * (w / wp) * (delta / wp) - (w / wp).^2);
end

