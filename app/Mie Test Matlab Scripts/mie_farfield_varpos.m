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
sf_layer = 2;
projector_padding = ceil((tf_layer + 1) / 2);

fs = c0 / 500e-9;
delay = 1 / fs;
ft = c0 / 500e-9;

rho = linspace(1000, 2000, 20);
phi = linspace(0, 2 * pi, 20);
th = linspace(0,  pi, 20);

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);
Ft = Ft(:);
Th = Th(:);
Phi = Phi(:);
Rho = Rho(:);
probes = [Th, Phi, Rho, Ft];

% far field probes
filename_probes = 'probes.in';
fileID = fopen(filename_probes, 'w');
fprintf(fileID, '%e %e %e\n', [-1, -1, -1] * dx);
fprintf(fileID, '%e %e %e\n', (dim + 1) * dx);
fprintf(fileID, '%d\n', size(probes, 1));
fprintf(fileID, '%e %e %e %e\n', probes');
fclose(fileID);

Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @Au;
er = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 60e-9;
size_param = K * a;

figure(1)

ft_samples = linspace(c0 / 1000e-9, c0 / 200e-9, 100);
plot(c0 ./ ft_samples / 1e-9, real(er_func(ft_samples * 2 * pi)), 'r-'), hold on
plot(c0 ./ ft_samples / 1e-9, -imag(er_func(ft_samples * 2 * pi)), 'b-')
xlabel('Wavelength [nm]')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, 'basic {\n');
fprintf(fileID, '%e %e\n', dt, dx);
fprintf(fileID, '%d %d %d\n', dim);
fprintf(fileID, '%d\n', sf_layer);
fprintf(fileID, '%d\n', tf_layer);
fprintf(fileID, '%d\n', PML_d);
fprintf(fileID, '%d\n', step);
fprintf(fileID, '%e %e\n', er_bg, ur_bg);
fprintf(fileID, '}\n');

% medium configuration
fprintf(fileID, 'medium 2 {\n');
% medium 0, background medium
fprintf(fileID, '{ ');
fprintf(fileID, '%e %e %e %e 0', er_bg, 0, ur_bg, 0);
fprintf(fileID, ' }\n');

% medium 1, scatterer medium
Au_print(fileID);
fprintf(fileID, '}\n');

% % geometry configuration
fprintf(fileID, 'geometry 1 {\n');
% geometry 0 gold sphere 10nm radius
fprintf(fileID, '{ ');
fprintf(fileID, 'sphere 1 %e %e %e %e', a, dim * dx / 2);
fprintf(fileID, ' }\n');
fprintf(fileID, '}\n');

% plane wave source
fprintf(fileID, 'source 1 {\n');
fprintf(fileID, '{ ');
fprintf(fileID, 'plane %d %d %s %e %e 0', projector_padding, dim(3) + projector_padding, 'r', fs, delay);
fprintf(fileID, ' }\n');
fprintf(fileID, '}\n');

% farfield probes
fprintf(fileID, 'farfield %s %s\n', filename_probes, 'output.out');
fclose(fileID);

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

