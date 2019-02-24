%% mie theory test, assuming vaccum background 
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

er_bg = 1;
ur_bg = 1;
PML_d = 6;
Sc = 0.5;
dx = 2e-9;
dt = Sc * dx / c0;
dim = [30, 30, 30];
step = 10000;
tf_layer = 2;
sf_layer = 3;
projector_padding = ceil((tf_layer + 1) / 2);
center = dim * dx / 2;
p2 = dim * dx;

fs = c0 / 500e-9;
delay = 1 / fs;
ft = c0 / 800e-9;

rho = linspace(1000, 2000, 20);
phi = linspace(0, 2 * pi, 20);
th = linspace(0,  pi, 20);

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);
Ft = Ft(:);
Th = Th(:);
Phi = Phi(:);
Rho = Rho(:);

ff_p1 = [-2, -2, -2] * dx;
ff_p2 = (dim + 2) * dx;

Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @Au;
[~, er] = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 30e-9;
size_param = K * a;

inhom.center = center;
inhom.len = [60e-9, 60e-9, 60e-9];
inhom.p1 = inhom.center - inhom.len / 2;
inhom.p2 = inhom.center + inhom.len / 2;
inhom.medium_idx1 = 1;
inhom.medium_idx2 = 0;
[inhom.x, inhom.y, inhom.z, inhom.dim, inhom.dx] = make_interp_region(inhom.p1, inhom.p2, 'dx', dx);
inhom.rho = (inhom.x - center(1)).^2 + (inhom.y - center(2)).^2 + (inhom.z - center(3)).^2 <= a^2;
inhom.rho = inhom.rho * 0.5;
inhom.type = 'inhom';
inhom.position = inhom.p1;
inhom.filename = 'inhom.in';

figure(1)

ft_samples = linspace(c0 / 1000e-9, c0 / 200e-9, 100);
[~, er_samples] = er_func(ft_samples * 2 * pi);
plot(c0 ./ ft_samples / 1e-9, real(er_samples), 'r-'), hold on
plot(c0 ./ ft_samples / 1e-9, -imag(er_samples), 'b-'), hold off
xlabel('Wavelength [nm]')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

disp('config.in created');
%%
basic.er_bg = er_bg;
basic.ur_bg = ur_bg;
basic.PML_d = PML_d;
basic.dx = dx;
basic.dt = dt;
basic.dim = dim;
basic.step = step;
basic.tf_layer = tf_layer;
basic.sf_layer = sf_layer;

medium{1} = Air();
medium{2} = Au();

% geometry{1} = struct('type', 'sphere', 'medium_idx', 1, 'radius', a, 'position', dim * dx / 2);
geometry{1} = inhom;

source{1} = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos', ...
    dim(3) + projector_padding, 'func_type', 'r', 'fp', fs, 'delay', delay, 'ref_pos', dim(3) / 2 * dx);

ff.lower_position = ff_p1;
ff.upper_position = ff_p2;
ff.theta = Th;
ff.phi = Phi;
ff.rho = Rho;
ff.freq = Ft;
ff.input_file = 'ff.in';
ff.output_file = 'output.out';

gen_config(basic, medium, geometry, source, 'farfield', ff, 'step_output', 1);

disp('config.in created');


%% numerical fields
call_exe('std_config')
data = load(ff.output_file);
make_complex = @(x, y) x + 1j * y;
ref_signal = load('reference.out');
ref_signal = reshape(ref_signal, 1, []);
t = (0:numel(ref_signal) - 1) * dt;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

Eth = make_complex(data(:, 1), data(:, 2)) ./ ref;
Ephi = make_complex(data(:, 3), data(:, 4)) ./ ref;
Hth = make_complex(data(:, 5), data(:, 6)) ./ ref;
Hphi = make_complex(data(:, 7), data(:, 8)) ./ ref;

disp('Numerical Fields Extracted');

%% theoretical fields
Eth_phy = zeros(size(m));
Ephi_phy = zeros(size(m));

for i = 1 : size(m, 1)
    res = Mie_S12(m(i), size_param(i), cos(Th(i)));
    S1 = res(1); 
    S2 = res(2);
    Eth_phy(i) = exp(1j * K(i) * Rho(i)) / (-1j * K(i) * Rho(i)) * cos(Phi(i)) * S2;
    Ephi_phy(i) = exp(1j * K(i) * Rho(i)) / (1j * K(i) * Rho(i)) * sin(Phi(i)) * S1;
end


%% corellation plots
figure(3)
plot(abs(Eth(:)), abs(Eth_phy(:)), '.'), hold on
plot(abs(Eth(:)), abs(Eth(:)), '-'), hold off
title('$|E_\theta|$', 'interpreter', 'latex','fontsize', 15)
axis equal
axis tight

figure(4)
plot(abs(Ephi(:)), abs(Ephi_phy(:)), '.'), hold on
plot(abs(Ephi(:)), abs(Ephi(:)), '-'), hold off
title('$|E_\phi|$', 'interpreter', 'latex','fontsize', 15)
axis equal
axis tight

