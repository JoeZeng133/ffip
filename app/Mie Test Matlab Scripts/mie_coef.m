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
tf_layer = 4;
sf_layer = 1;
projector_padding = ceil((tf_layer + 1) / 2);
a = 60e-9;
lam_min = 200e-9;
lam_max = 1000e-9;

lam = linspace(lam_min, lam_max, 50);
ft = c0 ./ lam;
Omega = 2 * pi * ft;
K = Omega / c0;

fs = c0 / (500e-9);
delay = 1/fs;
er_func = @Au;
[~, er] = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
size_param = K * a;

figure(1)
ft_samples = ft;
[~, er_samples] = er_func(ft_samples * 2 * pi);
plot(c0 ./ ft_samples / 1e-9, real(er_samples), 'r-'), hold on
plot(c0 ./ ft_samples / 1e-9, -imag(er_samples), 'b-')
xlabel('Wavelength [nm]')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

flux_p1 = [-1, -1, -1] * dx;
flux_p2 = (dim + 1) * dx;

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

geometry{1} = struct('type', 'sphere', 'medium_idx', 1, 'radius', a, 'position', dim * dx / 2);

source{1} = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos', ...
    dim(3) + projector_padding, 'func_type', 'r', 'fp', fs, 'delay', delay, 'ref_pos', dim(3) / 2 * dx);

flux.lower_position = flux_p1;
flux.upper_position = flux_p2;
flux.freq = ft;
flux.input_file = 'flux.in';
flux.output_file = 'flux.out';

gen_config(basic, medium, geometry, source, 'flux', flux);
disp('config.in created');

%% theoretical fields
c_scat_phy = zeros(size(m));
c_abs_phy = zeros(size(m));

for i = 1 : size(m, 1)
    res = Mie(m(i), size_param(i));
    c_scat_phy(i) = res(2);
    c_abs_phy(i) = res(3);
end

%% numerical fields
ref_signal = load('reference.out');
ref_signal = reshape(ref_signal, 1, []);
t = (0:numel(ref_signal) - 1) * dt;
ref = sum(ref_signal .* exp(-1j * Omega(:) * t), 2);
c_scat = load(flux.output_file);
c_scat = -c_scat ./ (0.5 * abs(ref).^2 / eta0 * pi * a^2);
disp('Numerical Fields Extracted');

%% comparison
figure(3)
plot(lam / 1e-9, c_scat, 'b*'), hold on
plot(lam / 1e-9, c_abs_phy, 'r-'), hold off
xlabel('Wavelength [nm]')
ylabel('Q_{abs}')
legend({'Simulation', 'Mie Theory'});
