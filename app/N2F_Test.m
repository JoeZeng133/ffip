%% far-field dipole test
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [10, 10, 10];
step = 600;

Np = 20;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, fp, 1/fp);

rho = linspace(1000, 2000, 11) * (Np * dx);
phi = linspace(0, 2 * pi, 10);
th = linspace(0 * pi,  pi, 10);
ft = linspace(0.5 * fp, 1.5 * fp, 20);

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);
pos = [Ft(:), Th(:), Phi(:), Rho(:)];

fileID = fopen('request.in', 'w');
fprintf(fileID, '%d\n', size(pos, 1));
fprintf(fileID, '%e %e %e %e\n', pos');
fclose(fileID);


%% theoretical fields
Omega = 2 * pi * Ft;
er =  1;
e = e0 * er;
ur = 1;
u = u0;

c = 1 ./ sqrt(e * u);
eta = sqrt(u ./ e);
K = Omega ./ c;

Hphi_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * K ./ Rho + 1 ./ Rho.^2) .* sin(Th);
Er_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (2 * eta ./ Rho.^2 + 2 ./ (1j * Omega .* e .* Rho.^3)) .* cos(Th);
Eth_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * Omega .* u ./ Rho + 1 ./ (1j * Omega .* e .* Rho.^3) + eta ./ Rho.^2) .* sin(Th);


%% simulated fields
data = load('data.out');
make_complex = @(x, y) x + 1j * y;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

Eth = make_complex(data(:, 1), data(:, 2)) ./ ref;
Ephi = make_complex(data(:, 3), data(:, 4)) ./ ref;
Hth = make_complex(data(:, 5), data(:, 6)) ./ ref;
Hphi = make_complex(data(:, 7), data(:, 8)) ./ ref;

%% comparisons
fprintf('Ephi over Eth = %e \n', max(abs(Ephi)) / max(abs(Eth)))
fprintf('Hth over Hphi = %e \n', max(abs(Hth)) / max(abs(Hphi)))

%% correlation comparisons
figure(3)
plot(abs(Eth(:)), abs(Eth_p(:)), '.')
xlabel('Numerical |E_\theta|', 'FontSize', 15)
ylabel('Theoretical |E_\theta|', 'FontSize', 15)
axis equal

figure(4)
plot(abs(Hphi(:)), abs(Hphi_p(:)), '.')
xlabel('Numerical |H_\phi|', 'FontSize', 15)
ylabel('Theoretical |H_\phi|', 'FontSize', 15)
axis equal
