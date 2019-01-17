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
er = er_func(Omega);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
size_param = K * a;

figure(1)
ft_samples = ft;
plot(c0 ./ ft_samples / 1e-9, real(er_func(ft_samples * 2 * pi)), 'r-'), hold on
plot(c0 ./ ft_samples / 1e-9, -imag(er_func(ft_samples * 2 * pi)), 'b-')
xlabel('Wavelength [nm]')
legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
axis tight

% scattering coefficients setup
flux_input_file = 'flux.in';
flux_output_file = 'flux.out';
fileID = fopen(flux_input_file, 'w');
fprintf(fileID, "%e %e %e\n", [-1, -1, -1] * dx);
fprintf(fileID, "%e %e %e\n", (dim + 1) * dx);
fprintf(fileID, '%d\n', numel(ft));
fprintf(fileID, '%e\n', ft);
fclose(fileID);

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
fprintf(fileID, "medium 1 {\n");
% medium Au
Au_print(fileID);
fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 1 {\n");
% geometry 0 gold sphere
fprintf(fileID, "{ ");
fprintf(fileID, "sphere 0 %e %e %e %e", a, dim * dx / 2);
fprintf(fileID, " }\n");
% geometry end
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, 'plane %d %d %s %e %e 0', projector_padding, dim(3) + projector_padding, 'r', fs, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% scattering coefficients calculated from poynting flux
fprintf(fileID, "flux %s %s\n", flux_input_file, flux_output_file);

disp('config.in created');

% theoretical fields
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
c_scat = load(flux_output_file);
c_scat = -c_scat ./ (0.5 * abs(ref).^2 / eta0 * pi * a^2);
disp('Numerical Fields Extracted');

%% comparison
figure(3)
plot(lam / 1e-9, c_scat, 'b*'), hold on
plot(lam / 1e-9, c_abs_phy, 'r-'), hold off
xlabel('Wavelength [nm]')
ylabel('Q_{abs}')
legend({'Simulation', 'Mie Theory'});
