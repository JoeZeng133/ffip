%% near-field dipole testclear
clc
clear
close all


e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 0.5;
dx = 1e-9;
dt = Sc * dx / c0;
step = 15000;
PML_d = 6;
dim = [22, 22, 22];
er_bg = 1;
ur_bg = 1;
tf_layer = 0;
sf_layer = 0;
projector_padding = ceil((tf_layer + 1) / 2);
center = dim * dx / 2;

lam_min = 200e-9;
lam_max = 1000e-9;
lam = linspace(lam_min, lam_max, 100);
fs = c0 / (500e-9);
delay = 1 / fs;
G = 1;

rho = 10 * dx;
phi = pi / 4;
th = pi / 4;
ft = c0 ./ lam;

[Rho, Phi, Th, Ft] = ndgrid(rho, phi, th, ft);
[Xo, Yo, Zo] = sph2cart(Phi, pi / 2 - Th, Rho);

%% theoretical fields
Omega = 2 * pi * Ft;
er_func = @Au;

figure(1)
ft_samples = ft;
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
basic.tf_layer = tf_layer;
basic.sf_layer = sf_layer;

medium{1} = Au();

geometry{1} = struct('type', 'box', 'medium_idx', 0, 'lower_position', [-1 -1 -1], 'upper_position', [1 1 1]);

source{1} = struct('type', 'dipole', 'filename', 'dipoles.in', 'x', center(1), 'y', center(2), 'z', center(3),...
    'amp', abs(G), 'fp', fs, 'func_type', double('r'), 'delay', delay, 'ctype', 4);

nf.x = Xo + center(1);
nf.y = Yo + center(2);
nf.z = Zo + center(3);
nf.freq = Ft;
nf.input_file = 'nf.in';
nf.output_file = 'output.out';

gen_config(basic, medium, geometry, source, 'nearfield', nf);

%% simulated fields
call_exe('std_config')
data = load(nf.output_file);
make_complex = @(x, y) x + 1j * y;

t = (0:step) * dt;
ref_signal = ricker(t, fs, delay);
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

E = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

proj_r = [sin(Th(:)) .* cos(Phi(:)), sin(Th(:)) .* sin(Phi(:)), cos(Th(:))];
proj_th = [cos(Th(:)) .* cos(Phi(:)), cos(Th(:)) .* sin(Phi(:)), -sin(Th(:))];
proj_phi = [-sin(Phi(:)), cos(Phi(:)), zeros(size(Phi(:)))];

Er = sum(E .* proj_r, 2) ./ ref;
Eth = sum(E .* proj_th, 2) ./ ref;
Hphi = sum(H .* proj_phi, 2) ./ ref;

%% comparisons
figure(2)
subplot(1, 2, 1)
plot(real(Er(:)), 'r'), hold on
plot(real(Er_p(:)), 'b'), hold off
subplot(1, 2, 2)
plot(imag(Er(:)), 'r'), hold on
plot(imag(Er_p(:)), 'b'), hold off

figure(3)
subplot(1, 2, 1)
plot(real(Eth(:)), 'r'), hold on
plot(real(Eth_p(:)), 'b'), hold off
subplot(1, 2, 2)
plot(imag(Eth(:)), 'r'), hold on
plot(imag(Eth_p(:)), 'b'), hold off

figure(4)
subplot(1, 2, 1)
plot(real(Hphi(:)), 'r'), hold on
plot(real(Hphi_p(:)), 'b'), hold off
subplot(1, 2, 2)
plot(imag(Hphi(:)), 'r'), hold on
plot(imag(Hphi_p(:)), 'b'), hold off