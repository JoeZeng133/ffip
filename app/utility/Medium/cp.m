function [res] = cp(w, A, phi, Omega, Gamma)
%cp pole, w in radius
%   Detailed explanation goes here
Omega = Omega * 2 * pi;
Gamma = Gamma * 2 * pi;
res = A * Omega * (exp(1j * phi) ./ (Omega + w - 1j * Gamma) + exp(-1j * phi) ./ (Omega - w + 1j * Gamma));
end

