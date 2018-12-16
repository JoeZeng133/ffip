clear
clc
close all

%% configuration input files
% constants
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

% basic configuration
er_bg = 1;
ur_bg = 1;
Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [50, 50, 50];
step = 600;
PMl_d = 6;

% excitation reference and frequency of interests (1 frequency)
Np = 30;                           
fp = c0 / (Np * dx);
d = 1 / fp;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, fp, d);
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * fp * t));

% medium for the inhomogeneous region
Omega = 2 * pi * fp;
lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
deybe = @(w, rel_e, tau) rel_e ./ (1 + 1j * w * tau);

er_const = 2 - 0.2i;
sig_const = -imag(er_const) * (2 * pi * fp * e0);
er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));

% Geometry for the inhomogeneous region
sphere_func = @(X, Y, Z, r, center) ((X - center(1)).^2 + (Y - center(2)).^2 + (Z - center(3)).^2) <= r^2;
geo_func = @(X, Y, Z) sphere_func(X, Y, Z, 10, dim / 2);

inhom_p1 = [15, 15, 15];
inhom_dim = [20, 20, 20];
inhom_p2 = inhom_p1 + inhom_dim;
inhom_x = inhom_p1(1):inhom_p2(1);
inhom_y = inhom_p1(2):inhom_p2(2);
inhom_z = inhom_p1(3):inhom_p2(3);
[inhom_X, inhom_Y, inhom_Z] = ndgrid(inhom_x, inhom_y, inhom_z);
inhom_pos = [inhom_X(:), inhom_Y(:), inhom_Z(:)];
num_inhom_pos = size(inhom_pos, 1);
rho_target = geo_func(inhom_X, inhom_Y, inhom_Z) * 1;
rho = ones(size(rho_target));                          %start with rho=1
rho = rho(:);

probes_pos = [15 15 40];                                %point to maximize
num_probes = 1;
%% simulated fields
make_complex = @(x, y) x + 1j * y;

E_target_probes = [0, 0, 0];
H_target_probes = [0, 0, 0];

save('case_configuration');
disp('objective fields saved');
