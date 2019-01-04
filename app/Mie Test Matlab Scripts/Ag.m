function res = Ag(w)
drude = @(w, wd, gamma) wd^2 ./ (1j * w * gamma - w.^2);
cp = @(w, A, phi, Omega, Gamma) A * Omega * (exp(1j * phi) ./ (Omega + w - 1j * Gamma) + exp(-1j * phi) ./ (Omega - w + 1j * Gamma));
res = 1.4447+ drude(w, 1.3280e16, 9.1269e13) + cp(w, -1.5951, 3.1288, 8.2749e15, 5.177e15) + cp(w, 0.25261, -1.5066, 6.1998e15, 5.4126e14);
end