%% mie theory test, assuming vaccum background 
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 1 / sqrt(3);
dt = 2e-17 / 50;
dx = c0 * dt / Sc;
dim = [50, 50, 50];
step = 1500;

Np = 30;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, fp, 1/fp);

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
ouput = [Ft, Th, Phi, Rho];

fileID = fopen('../request.in', 'w');
fprintf(fileID, '%d\n', size(ouput, 1));
fprintf(fileID, '%e %e %e %e\n', ouput');
fclose(fileID);

%%
Omega = 2 * pi * Ft;
K = Omega / c0;

lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
er_func = @(w) (1 + lorentz(w, 0.8, 4e16, 1e16) + lorentz(w, 0.5, 6e16, 1e16));
er = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 10 * dx;
size_param = K * a;

% inhomogeneous geometry file generation

fileID = fopen('../geometry.in', 'w');
vspan = 30;
center = dim * dx / 2;
st = center - 0.5 * vspan * dx;
v = (0:vspan) * dx + st(1);
[X, Y, Z] = ndgrid(v, v, v);
V = ((X - center(1)).^2 + (Y - center(2)).^2 + (Z - center(3)).^2) < a^2;
fprintf(fileID, '%d %d %d\n', [vspan+1, vspan+1, vspan+1]);
fprintf(fileID, '%e %e\n', [st(1), dx]);
fprintf(fileID, '%e %e\n', [st(2), dx]);
fprintf(fileID, '%e %e\n', [st(2), dx]);
fprintf(fileID, '%e ', V(:));
fclose(fileID);

figure(1)
ft_samples = linspace(0.5 * fp, 1.5 * fp, 100);
plot(ft_samples / 1e15, real(er_func(ft_samples * 2 * pi)), 'r-'), hold on
plot(ft_samples / 1e15, -imag(er_func(ft_samples * 2 * pi)), 'b-')
xlabel('Frequency (PHz)')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight


%% numerical fields
data = load('../data.out');
make_complex = @(x, y) x + 1j * y;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

Eth = make_complex(data(:, 1), data(:, 2)) ./ ref;
Ephi = make_complex(data(:, 3), data(:, 4)) ./ ref;
Hth = make_complex(data(:, 5), data(:, 6)) ./ ref;
Hphi = make_complex(data(:, 7), data(:, 8)) ./ ref;

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

