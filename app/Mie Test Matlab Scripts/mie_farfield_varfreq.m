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
dim = [60, 60, 60];
step = 10000;
tf_layer = 2;
sf_layer = 3;
projector_padding = ceil((tf_layer + 1) / 2);

fs = c0 / 500e-9;
delay = 1 / fs;

rho = 1000;
phi = pi / 4;
th = pi / 4;
ft = linspace(c0 / 200e-9, c0 / 1000e-9, 100);

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);
Ft = Ft(:);
Th = Th(:);
Phi = Phi(:);
Rho = Rho(:);

% far field probes
ff_p1 = [-2, -2, -2] * dx;
ff_p2 = (dim + 2) * dx;

Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @Au;
[~, er] = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 60e-9;
size_param = K * a;

figure(1)

ft_samples = ft;
[~, er_samples] = er_func(ft_samples * 2 * pi);
plot(c0 ./ ft_samples / 1e-9, real(er_samples), 'r-'), hold on
plot(c0 ./ ft_samples / 1e-9, -imag(er_samples), 'b-')
xlabel('Wavelength [nm]')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

%% writing configuration
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

geometry{1} = struct('type', 'sphere', 'medium_idx', 1, 'radius', a, 'position', dim * dx / 2);

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

gen_config(basic, medium, geometry, source, 'farfield', ff);

disp('config.in created');
%% numerical fields
data = load('output.out');
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


%% plots against one parameters (frequency)

figure(3)
plot(c0 ./ ft / 1e-9, abs(Eth(:)), 'sq'), hold on
plot(c0 ./ ft / 1e-9, abs(Eth_phy(:)))
xlabel('Wavelength [nm]')
legend({'$Numerical |E_\theta|$', '$Analytical |E_\theta|$'}, 'interpreter', 'latex','fontsize', 15)

figure(4)
plot(c0 ./ ft / 1e-9, abs(Ephi(:)), 'sq-'), hold on
plot(c0 ./ ft / 1e-9, abs(Ephi_phy(:)))
xlabel('Wavelength [nm]')
legend({'$Numerical |E_\phi|$', '$Analytical |E_\phi|$'}, 'interpreter', 'latex','fontsize', 15)

