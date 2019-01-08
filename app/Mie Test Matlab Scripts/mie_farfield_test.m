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
Sc = 1 / sqrt(3);
dt = 2e-17 / 50;
dx = c0 * dt / Sc;
dim = [25, 25, 25];
step = 600;

Np = 30;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
d = 1/fp;
ref_signal = ricker(t, fp, d);

% rho = linspace(1000, 2000, 11) * (Np * dx);
% phi = linspace(0, 2 * pi, 10);
% th = linspace(0 * pi,  pi, 10);
rho = 1;
phi = pi / 4;
th = pi / 4;
ft = linspace(0.5 * fp, 1.5 * fp, 100);

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);
Ft = Ft(:);
Th = Th(:);
Phi = Phi(:);
Rho = Rho(:);
probes = [Th, Phi, Rho, Ft];

% write to probes file
filename_probes = 'probes.in';
fileID = fopen(filename_probes, "w");
fprintf(fileID, "%e %e %e\n", [-1, -1, -1] * dx);
fprintf(fileID, "%e %e %e\n", (dim + 1) * dx);
fprintf(fileID, "%d\n", size(probes, 1));
fprintf(fileID, "%e %e %e %e\n", probes');
fclose(fileID);

lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);


Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @(w) (1 + lorentz(w, 0.8, 4e16, 1e16) + lorentz(w, 0.5, 6e16, 1e16));
% er_func = @(w) (1 + drude(w, 2e14, 0.5e14) + lorentz(w, 3, 5e14, 1e14));
er = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 10 * dx;
size_param = K * a;

figure(1)
ft_samples = linspace(0.5 * fp, 1.5 * fp, 100);
plot(ft_samples / 1e15, real(er_func(ft_samples * 2 * pi)), 'r-'), hold on
plot(ft_samples / 1e15, -imag(er_func(ft_samples * 2 * pi)), 'b-')
xlabel('Frequency (PHz)')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

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
fprintf(fileID, "medium 2 {\n");
% medium 0, background medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
fprintf(fileID, " }\n");

% medium 1, scatterer medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 2\n", er_bg, 0, ur_bg, 0);
fprintf(fileID, "{ Lorentz %e %e %e }\n", 0.8, 4e16, 1e16);
fprintf(fileID, "{ Lorentz %e %e %e }\n", 0.5, 6e16, 1e16);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% % geometry configuration
fprintf(fileID, "geometry 1 {\n");
% geometry 0 gold sphere 60nm
fprintf(fileID, "{ ");
fprintf(fileID, "sphere 1 %e %e %e %e", a, dim * dx / 2);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "plane %d %e %e", dim(3), fp, d);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% farfield probes
fprintf(fileID, "farfield %s %s\n", filename_probes, 'output.out');
fclose(fileID);

disp('config.in created');
%% numerical fields
data = load('output.out');
make_complex = @(x, y) x + 1j * y;
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
plot(Ft(:)/1e15, abs(Eth(:)), 'sq'), hold on
plot(Ft(:)/1e15, abs(Eth_phy(:)))
xlabel('Frequency (PHz)')
legend({'$Numerical |E_\theta|$', '$Analytical |E_\theta|$'}, 'interpreter', 'latex','fontsize', 15)


figure(4)
plot(Ft(:)/1e15, abs(Ephi(:)), 'sq-'), hold on
plot(Ft(:)/1e15, abs(Ephi_phy(:)))
xlabel('Frequency (PHz)')
legend({'$Numerical |E_\phi|$', '$Analytical |E_\phi|$'}, 'interpreter', 'latex','fontsize', 15)


%% corellation plots
% figure(3)
% plot(abs(Eth(:)), abs(Eth_phy(:)), 'sq')
% title('$|E_\theta|$', 'interpreter', 'latex','fontsize', 15)
% axis equal
% axis tight
% 
% figure(4)
% plot(abs(Ephi(:)), abs(Ephi_phy(:)), 'sq')
% title('$|E_\phi|$', 'interpreter', 'latex','fontsize', 15)
% axis equal
% axis tight

