close all
clear
clc

load("adjoint_results.mat");

%%
trap_weights = @(dim) [0.5, ones([1, dim - 2]), 0.5];
[w_X, w_Y, w_Z] = ndgrid(trap_weights(inhom_dim(1) + 1), trap_weights(inhom_dim(2) + 1), trap_weights(inhom_dim(3) + 1));
w = w_X .* w_Y .* w_Z;
w = w(:) * dx^3;
Se = 1j * (2 * pi * fp) * sum(E_forward_inhom .* E_adjoint_inhom, 2);
A = real(Se .* w * e0 * (er_const - er_bg));
A = reshape(A, inhom_dim + 1);
target = sum(sum((abs(E_forward_probes) - abs(E_target_probes)).^2, 2), 1);

%%
drho = linprog(A(:), [], [], [], [], -rho(:), 1 - rho(:));
