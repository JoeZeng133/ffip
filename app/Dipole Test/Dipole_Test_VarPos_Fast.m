%% near-field dipole testclear
clc
clear
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 0.5;
dt = 7.2e-17;
dx = c0 * dt / Sc;
step = 1000;
PML_d = 6;
dim = [22, 22, 22];
er_bg = 1;
ur_bg = 1;
tf_layer = 0;
sf_layer = 0;
projector_padding = ceil((tf_layer + 1) / 2);
center = dim * dx / 2;

l = 30 * dx;
ft = c0 / l;
fs = c0 / l;
delay = 1 / fs;
G1 = 1 - 0.3j;
G2 = 2 + 1j;
G = G1 + G2;

rho = linspace(5, 10, 10) * dx;
phi = linspace(0, 2 * pi, 10);
th = linspace(0, pi, 10);

[Rho, Phi, Th, Ft] = ndgrid(rho, phi, th, ft);
[Xo, Yo, Zo] = sph2cart(Phi, pi / 2 - Th, Rho);

% wriet to dipole source files
dipoles.type = 'dipole';
dipoles.amp = abs([G1 G2]);
dipoles.ctype = [4 4];
dipoles.fp = [fs fs];
dipoles.x = [center(1) center(1)];
dipoles.y = [center(2) center(2)];
dipoles.z = [center(3) center(3)];
dipoles.func_type = double('r') * [1 1];
dipoles.delay = delay - angle([G1 G2]) / (2 * pi * ft);
dipoles.filename = 'dipoles.in';

%% theoretical fields
Omega = 2 * pi * Ft;
lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
deybe = @(w, rel_e, tau) rel_e ./ (1 + 1j * w * tau);

er_func = @fic1;

figure(1)
ft_samples = c0 ./ linspace(100e-9, 2000e-9, 100);
[~, er_samples] = er_func(2 * pi * ft_samples);
plot(3e8 ./ ft_samples / 1e-9, real(er_samples)), hold on
plot(3e8 ./ ft_samples / 1e-9, -imag(er_samples), '--'), hold off
legend({['Re ', char(949), '_r'], ['Im ', char(949), '_r']}, 'FontSize', 15)
xlabel('\lambda [nm]')

[~, er] = er_func(Omega);
e = e0 * er;

u = u0;
c = 1 ./ sqrt(e * u);
eta = sqrt(u ./ e);
K = Omega ./ c;


Hphi_p = G * 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * K ./ Rho + 1 ./ Rho.^2) .* sin(Th);
Er_p = G * 1 / (4 * pi) * exp(-1j * K .* Rho) .* (2 * eta ./ Rho.^2 + 2 ./ (1j * Omega .* e .* Rho.^3)) .* cos(Th);
Eth_p = G * 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * Omega .* u ./ Rho + 1 ./ (1j * Omega .* e .* Rho.^3) + eta ./ Rho.^2) .* sin(Th);

%% write configurations
basic.er_bg = er_bg;
basic.ur_bg = ur_bg;
basic.PML_d = PML_d;
basic.dx = dx;
basic.dt = dt;
basic.dim = dim;
basic.step = step;
basic.tf_layer = 1;
basic.sf_layer = 0;

medium{1} = fic1();

geometry{1} = struct('type', 'box', 'medium_idx', 0, 'lower_position', [-1 -1 -1], 'upper_position', [1 1 1]);

source{1} = dipoles;

nf.x = Xo + center(1);
nf.y = Yo + center(2);
nf.z = Zo + center(3);
nf.freq = Ft;
nf.input_file = 'nf.in';
nf.output_file = 'output.out';

gen_config(basic, medium, geometry, source, 'nearfield', nf, 'num_proc', 2, 'step_output', 1);

%% simulated fields
call_exe('std_config')
data = load('output.out');
make_complex = @(x, y) x + 1j * y;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);

t = (0:step) * dt;
ref_signal = ricker(t, fs, delay);
% ref_signal = sin(2 * pi * fp * (t - delay));
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

E = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

proj_r = [sin(Th(:)) .* cos(Phi(:)), sin(Th(:)) .* sin(Phi(:)), cos(Th(:))];
proj_th = [cos(Th(:)) .* cos(Phi(:)), cos(Th(:)) .* sin(Phi(:)), -sin(Th(:))];
proj_phi = [-sin(Phi(:)), cos(Phi(:)), zeros(size(Phi(:)))];

Er = sum(E .* proj_r, 2) ./ ref;
Eth = sum(E .* proj_th, 2) ./ ref;
Hphi = sum(H .* proj_phi, 2) ./ ref;

%% correlation comparisons
figure
subplot(1, 2, 1)
plot(real(Er(:)), real(Er_p(:)), '.'), hold on
plot(real(Er(:)), real(Er(:))), hold off
axis equal
axis tight
title('Re E_r')

subplot(1, 2, 2)
plot(imag(Er(:)), imag(Er_p(:)), '.'), hold on
plot(imag(Er(:)), imag(Er(:))), hold off
axis equal
axis tight
title('Im E_r')

figure
subplot(1, 2, 1)
plot(real(Eth(:)), real(Eth_p(:)), '.'), hold on
plot(real(Eth(:)), real(Eth(:))), hold off
axis equal
axis tight
title('Re E_\theta')

subplot(1, 2, 2)
plot(imag(Eth(:)), imag(Eth_p(:)), '.'), hold on
plot(imag(Eth(:)), imag(Eth(:))), hold off
axis equal
axis tight
title('Im E_\theta')

figure
subplot(1, 2, 1)
plot(real(Hphi(:)), real(Hphi_p(:)), '.'), hold on
plot(real(Hphi(:)), real(Hphi(:))), hold off
axis equal
axis tight
title('Re H_\phi')

subplot(1, 2, 2)
plot(imag(Hphi(:)), imag(Hphi_p(:)), '.'), hold on
plot(imag(Hphi(:)), imag(Hphi(:))), hold on
axis equal
axis tight
title('Im H_\phi')




