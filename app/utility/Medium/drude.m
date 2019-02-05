function [res] = drude(w, fd, gamma)
%drude pole, fd gamma in Hz
%   Detailed explanation goes here
wd = 2 * pi * fd;
gamma = 2 * pi * gamma;
res = wd^2 ./ (1j * w * gamma - w.^2);
end

