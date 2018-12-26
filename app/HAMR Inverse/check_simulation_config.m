%% check simulation setup using a Mie theory test
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

er_bg = 1;
ur_bg = 1;
c = c0 / sqrt(er_bg * ur_bg);

PML_d = 6;
dt = 4.5e-18;           
dx = 2.5e-9;            %2.5nm discretization
Sc = c * dt / dx;       %Sc < 1/sqrt(3) = 0.5774
dim = [50, 50, 50];
step = 1000;

fp = 5.4e14;            %540 THz 
Np = c / (fp * dx);     %Wavelength in background medium [normalized to dx]
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
delay = 1 / fp;
ref_signal = ricker(t, fp, delay);

% rho = linspace(1000, 2000, 11) * (Np * dx);
% phi = linspace(0, 2 * pi, 10);
% th = linspace(0 * pi,  pi, 10);
rho = 1;
phi = pi / 4;
th = pi / 4;
lam_min = 400e-9;
lam_max = 900e-9;
lam = linspace(lam_min, lam_max, 100);
ft = c ./ lam;

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);
Ft = Ft(:);
Th = Th(:);
Phi = Phi(:);
Rho = Rho(:);
probes = [Th, Phi, Rho, Ft];

% write to probes file
filename_probes = 'probes.in';
fileID = fopen(filename_probes, "w");
fprintf(fileID, "%d\n", size(probes, 1));
fprintf(fileID, "%e %e %e %e\n", probes');
fclose(fileID);

lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
 
Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @(w) (1 + drude(w, 1.323e16 / (2 * pi), 1.26e14 / (2 * pi)));   %gold for l > 750nm
er = er_func(Omega);

figure(1)
plot(lam / 1e-9, real(er), 'r-'), hold on
plot(lam / 1e-9, -imag(er), 'b-')
xlabel('\lambda [nm]')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

%%
% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 30e-9;
size_param = K * a;

% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "%d\n", PML_d);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 1 {\n");
% medium 0, Gold
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 1\n", er_bg, 0, ur_bg, 0);
fprintf(fileID, "{ Drude %e %e }\n", 1.323e16 / (2 * pi), 1.26e14 / (2 * pi));
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 1 {\n");
% geometry 0 gold sphere 60nm
fprintf(fileID, "{ ");
fprintf(fileID, "sphere 0 %e %e %e %e", a, dim * dx / 2);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "eigen %d %e %e", dim(3), fp, delay);
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
    res = Mie2_S12(m(i)^2, 1, size_param(i), cos(Th(i)));
    S1 = res(1); 
    S2 = res(2);
    Eth_phy(i) = exp(1j * K(i) * Rho(i)) / (-1j * K(i) * Rho(i)) * cos(Phi(i)) * S2;
    Ephi_phy(i) = exp(1j * K(i) * Rho(i)) / (1j * K(i) * Rho(i)) * sin(Phi(i)) * S1;
end


%% plots against one parameters (frequency)

figure(3)
plot(lam / 1e-9, abs(Eth(:)), 'sq'), hold on
plot(lam / 1e-9, abs(Eth_phy(:)))
xlabel('\lambda [nm]')
legend({'$Numerical |E_\theta|$', '$Analytical |E_\theta|$'}, 'interpreter', 'latex','fontsize', 15)


figure(4)
plot(lam / 1e-9, abs(Ephi(:)), 'sq-'), hold on
plot(lam / 1e-9, abs(Ephi_phy(:)))
xlabel('\lambda [nm]')
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

