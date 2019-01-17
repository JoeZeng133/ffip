function res = Au2(w)
drude = @(w, wd, gamma) wd^2 ./ (1j * w * gamma - w.^2);
lorentz = @(w, rel_e, wp, gamma) rel_e * wp^2 ./ ((wp)^2 - w.^2 + 1j * w * gamma);
% cp = @(w, A, phi, Omega, Gamma) A * Omega * (exp(1j * phi) ./ (Omega + w - 1j * Gamma) + exp(-1j * phi) ./ (Omega - w + 1j * Gamma));
res = 5.9673 + drude(w, 2113.6e12 * 2 * pi, 15.92e12 * 2 * pi) + lorentz(w, 1.09, 650.07e12 * 2 * pi, 104.86e12 * 2 * pi);
end