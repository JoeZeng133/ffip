clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
sim_dim = [30, 30, 30];
inhom_dim = [20, 20, 20];
p1 = 0.5 * (sim_dim - inhom_dim);
p2 = 0.5 * (sim_dim + inhom_dim);
num_samples = 200;
step = 600;

Np = 20;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
ricker = @(t, d, fp) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, 0, fp);

[X, Y, Z] = ndgrid(0:sim_dim(1), 0:sim_dim(2), 0:sim_dim(3));
X = X(X < p1(1) | X > p2(1));
Y = Y(Y < p1(2) | Y > p2(2));
Z = Z(Z < p1(3) | Z > p2(3));

output = datasample([[X(:), Y(:), Z(:)] * dx, fp * ones(size(X(:)))], num_samples, 1, 'Replace', false);
fileID = fopen('request.in', 'w');
fprintf(fileID, '%d\n', size(output, 1));
fprintf(fileID, '%e %e %e %e\n', output');
fclose(fileID);

fileID = fopen('geometry.in', 'w');
fprintf(fileID, '1\n');                         %1 geometry
fprintf(fileID, 'sphere\n');                     %inhomogeneous region

fclose(fileID);

%% theoretical fields
Omega = 2 * pi * Ft;
lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
deybe = @(w, rel_e, tau) rel_e ./ (1 + 1j * w * tau);

er_const = 2 - 0.2i;
sig_const = -imag(er_const) * (2 * pi * fp * e0);
er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));
% rel_e = (real(er_const)) * (1 + imag(er_const)^2 / real(er_const)^2);
% tau = -imag(er_const) / real(er_const) / (2 * pi * fp);
% er_func = @(w) (1 + deybe(w, rel_e, tau));
% er_func = @(w) (1 + drude(w, 5e14, 2e14));
% er_func = @(w) (1 + lorentz(w, 3, 3e14, 0.5e14) + lorentz(w, 3, 5e14, 1e14));
% er_func = @(w) (1 + drude(w, 2e14, 0.5e14) + lorentz(w, 3, 5e14, 1e14));

testf = linspace(0.5 * fp, 2 * fp, 1000);
tester = er_func(2 * pi * testf);

figure(1)
plot(testf / (1e12), real(tester)), hold on
plot(testf / (1e12), -imag(tester), '--'), hold off
legend({['Re ', char(949), '_r'], ['Im ', char(949), '_r']}, 'FontSize', 15)
xlabel('f [THz]')



%% simulated fields
data = load('data.out');
make_complex = @(x, y) x + 1j * y;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

E = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];


