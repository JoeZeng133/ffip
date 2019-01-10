function res = Au2(w)
drude = @(w, wd, gamma) wd^2 ./ (1j * w * gamma - w.^2);
cp = @(w, A, phi, Omega, Gamma) A * Omega * (exp(1j * phi) ./ (Omega + w - 1j * Gamma) + exp(-1j * phi) ./ (Omega - w + 1j * Gamma));
res = 7.1431 + drude(w, 1.3202e16, 1.0805e14);
end