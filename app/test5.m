%% near-field dipole test
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

er = 1;
ur = 1;
e = e0 * er;
u = u0 * ur;
c = 1 / sqrt(e * u);
eta = sqrt(u / e);


dt = 2e-17 / 50;
dx = c * dt / (1 / sqrt(3));
dim = [30, 30, 30];
step = 600;

Np = 20;                            %center frequency of the rickerwavelet
fp = c / (Np * dx);
ricker = @(t, d, fp) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, 0, fp);

rho = 12 * dx;
phi = pi / 4;
th = linspace(0, pi / 2, 20);
ft = linspace(0.8 * fp, 1.2 * fp, 20);

[Rho, Phi, Th, Ft] = ndgrid(rho, phi, th, ft);
[Xo, Yo, Zo] = sph2cart(Phi, pi / 2 - Th, Rho);

output = [[Xo(:), Yo(:), Zo(:)] + dim * dx / 2, Ft(:)];

fileID = fopen('request.in', 'w');
fprintf(fileID, '%d\n', size(output, 1));
fprintf(fileID, '%e %e %e %e\n', output');
fclose(fileID);


%% theoretical fields
Omega = 2 * pi * Ft;
K = Omega ./ c;

Hphi_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * K ./ Rho + 1 ./ Rho.^2) .* sin(Th);
Er_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (2 * eta ./ Rho.^2 + 2 ./ (1j * Omega * e .* Rho.^3)) .* cos(Th);
Eth_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * Omega .* u ./ Rho + 1 ./ (1j * Omega * e .* Rho.^3) + eta ./ Rho.^2) .* sin(Th);


%% simulated fields
data = load('data.out');
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
figure(1)
subplot(1, 2, 1)
plot(real(Er(:)), real(Er_p(:)), '*')
axis equal

subplot(1, 2, 2)
plot(imag(Er(:)), imag(Er_p(:)), '*')
axis equal
% plot(imag(Er(:)), '*'), hold on
% plot(imag(Er_p(:)))

figure(2)
subplot(1, 2, 1)
% plot(real(Eth(:)), '*'), hold on
% plot(real(Eth_p(:))), hold off
plot(real(Eth(:)), real(Eth_p(:)), '*')
axis equal

subplot(1, 2, 2)
% plot(imag(Eth(:)), '*'), hold on
% plot(imag(Eth_p(:))), hold off
plot(imag(Eth(:)), imag(Eth_p(:)), '*')
axis equal


