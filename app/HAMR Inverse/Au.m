function res = Au(w)
drude = @(w, wd, gamma) wd^2 ./ (1j * w * gamma - w.^2);
cp = @(w, A, phi, Omega, Gamma) A * Omega * (exp(1j * phi) ./ (Omega + w - 1j * Gamma) + exp(-1j * phi) ./ (Omega - w + 1j * Gamma));
res = 1.1431 + drude(w, 1.3202e16, 1.0805e14) + cp(w, 0.26698, -1.2371, 3.8711e15, 4.4642e14) + cp(w, 3.0834, -1.0968, 4.1684e15, 2.3555e15);
end