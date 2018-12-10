close all
clear
clc

load('case_configuration.mat');
load('adjoint_results.mat');
load('forward_results.mat');

%%
trap_weights = @(dim) [0.5, ones([1, dim - 2]), 0.5];
[w_X, w_Y, w_Z] = ndgrid(trap_weights(inhom_dim(1) + 1), trap_weights(inhom_dim(2) + 1), trap_weights(inhom_dim(3) + 1));
w = w_X .* w_Y .* w_Z;
w = w(:) * dx^3;
S_sigma = 1j * (2 * pi * fp) * sum(E_inhom_forward .* E_inhom_adjoint, 2);
A = real(S_sigma .* w * e0 * (er_const - er_bg));
A = reshape(A, inhom_dim + 1);

%%
surf(A(:, :, 10))