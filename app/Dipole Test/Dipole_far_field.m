%% near-field dipole testclear
clc
close all


e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [30, 30, 30];
step = 600;
PML_d = 6;
er_bg = 1;
ur_bg = 1;

Np = 20;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
delay = 1 / fp;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, fp, delay);
G = 1 + 2j;                         %the point source is G

rho = linspace(10, 15, 10) * dx;
phi = linspace(0, 2 * pi, 10);
th = linspace(1 / 4 * pi, 3 / 4 * pi, 10);
ft = fp;

[Rho, Phi, Th, Ft] = ndgrid(rho, phi, th, ft);
[Xo, Yo, Zo] = sph2cart(Phi, pi / 2 - Th, Rho);

probes = [[Xo(:), Yo(:), Zo(:)] + dim * dx / 2, Ft(:)];

% write to probes file
filename_probes = 'probes.in';
fileID = fopen(filename_probes, "w");
fprintf(fileID, "%d\n", size(probes, 1));
fprintf(fileID, "%e %e %e %e\n", probes');
fclose(fileID);

% wriet to dipole source files
filename_dipoles = 'dipoles.in';
fileID = fopen(filename_dipoles, 'w');
fprintf(fileID, '1\n');
fprintf(fileID, '%e %e %e %e %e %e 4', dim * dx / 2, abs(G), fp, delay - angle(G) / (2 * pi * fp));
fclose(fileID);

%% theoretical fields
Omega = 2 * pi * Ft;
lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
deybe = @(w, rel_e, tau) rel_e ./ (1 + 1j * w * tau);

% er_const = 2 - 0.2i;
% sig_const = -imag(er_const) * (2 * pi * fp * e0);
% er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));
% er_func = @(w) (1 + deybe(w, rel_e, tau));
er_func = @(w) (1 + drude(w, 2e14, 0.5e14));
% er_func = @(w) (1 + lorentz(w, 3, 3e14, 0.5e14) + lorentz(w, 3, 5e14, 1e14));
% er_func = @(w) (1 + drude(w, 2e14, 0.5e14) + lorentz(w, 3, 5e14, 1e14));
% er_func = @(w) (1 + lorentz(w, 3, 5e14, 1e14));

testf = linspace(0.5 * fp, 2 * fp, 1000);
tester = er_func(2 * pi * testf);

figure(1)
plot(testf / (1e12), real(tester)), hold on
plot(testf / (1e12), -imag(tester), '--'), hold off
legend({['Re ', char(949), '_r'], ['Im ', char(949), '_r']}, 'FontSize', 15)
xlabel('f [THz]')

er =  er_func(Omega);
e = e0 * er;

u = u0;
c = 1 ./ sqrt(e * u);
eta = sqrt(u ./ e);
K = Omega ./ c;


Hphi_p = G * 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * K ./ Rho + 1 ./ Rho.^2) .* sin(Th);
Er_p = G * 1 / (4 * pi) * exp(-1j * K .* Rho) .* (2 * eta ./ Rho.^2 + 2 ./ (1j * Omega .* e .* Rho.^3)) .* cos(Th);
Eth_p = G * 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * Omega .* u ./ Rho + 1 ./ (1j * Omega .* e .* Rho.^3) + eta ./ Rho.^2) .* sin(Th);

%% write configurations
% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 2);
fprintf(fileID, "%d\n", PML_d);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 1{\n");
% medium 0, scatterer medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 1\n", er_bg, 0, ur_bg, 0);
fprintf(fileID, "{ Drude %e %e }\n", 2e14, 0.5e14);
% fprintf(fileID, "{ Lorentz %e %e %e }\n", 3, 5e14, 1e14);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 1 {\n");
% geometry 0, the box region that covers all the fields
fprintf(fileID, "{ ");
fprintf(fileID, "box 0 -1 -1 -1 1 1 1");
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% dipole sources
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "dipole %s", filename_dipoles);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% nearfield probes
fprintf(fileID, "nearfield %s %s\n", filename_probes, 'output.out');
fclose(fileID);

disp('config.in created');
%% simulated fields
data = load('output.out');
make_complex = @(x, y) x + 1j * y;
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
% figure(2)
% plot(abs(Er(:)), 'r'), hold on
% plot(abs(Er_p(:)), 'b')

%% correlation comparisons
figure(3)
subplot(1, 2, 1)
plot(abs(Er(:)), abs(Er_p(:)), '*')
axis equal
axis tight
title('Re E_r')

subplot(1, 2, 2)
plot(imag(Er(:)), imag(Er_p(:)), '*')
axis equal
axis tight
title('Im E_r')

figure(4)
subplot(1, 2, 1)
plot(abs(Eth(:)), abs(Eth_p(:)), '*')
axis equal
axis tight
title('Re E_\theta')

subplot(1, 2, 2)
plot(imag(Eth(:)), imag(Eth_p(:)), '*')
axis equal
axis tight
title('Im E_\theta')

figure(5)
subplot(1, 2, 1)
plot(abs(Hphi(:)), abs(Hphi_p(:)), '*')
axis equal
axis tight
title('Re H_\phi')

subplot(1, 2, 2)
plot(imag(Hphi(:)), imag(Hphi_p(:)), '*')
axis equal
axis tight
title('Im H_\phi')


