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

lam_min = 200e-9;
lam_max = 1000e-9;
lam = linspace(lam_min, lam_max, 100);
fp = c0 / (500e-9);
delay = 1 / fp;
G = 1;

rho = 10 * dx;
phi = pi / 4;
th = pi / 4;
ft = c0 ./ lam;

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
fprintf(fileID, '%e %e %e %e %s %e %e 4\n', dim * dx / 2, abs(G), 'r', fp, delay);
fclose(fileID);

%% theoretical fields
Omega = 2 * pi * Ft;
lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
deybe = @(w, rel_e, tau) rel_e ./ (1 + 1j * w * tau);

er_func = @Au;
testf = ft;
tester = er_func(2 * pi * testf);

figure(1)
plot(3e8 ./ testf / 1e-9, real(tester)), hold on
plot(3e8 ./ testf / 1e-9, -imag(tester), '--'), hold off
legend({['Re ', char(949), '_r'], ['Im ', char(949), '_r']}, 'FontSize', 15)
xlabel('\lambda [nm]')

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
fprintf(fileID, "%d\n", 0);
fprintf(fileID, "%d\n", 0);
fprintf(fileID, "%d\n", PML_d);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 1{\n");
% Au
Au_print(fileID);
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
% call_exe('std_config')
data = load('output.out');
make_complex = @(x, y) x + 1j * y;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);

t = (0:step) * dt;
ref_signal = ricker(t, fp, delay);
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