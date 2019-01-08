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
dim = [32, 32, 32];
step = 3000;

lam_min = 200e-9;
lam_max = 1000e-9;
lam = linspace(lam_min, lam_max, 50);
ft = c0 ./ lam;
Omega = 2 * pi * ft;
K = Omega / c0;
er_func = @Au;
er = er_func(Omega);

fp = c0 / (500e-9);
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
delay = 1/fp;
ref_signal = ricker(t, fp, delay);
ref = sum(ref_signal .* exp(-1j * Omega(:) * t), 2);

% sphere parameters
m = sqrt(conj(er(:))); %exp(-jwt) dependence, use the conjugate
a = 30e-9;
size_param = K * a;

lam_plot = linspace(200e-9, 1000e-9, 100);
ft_plot = c0 ./ lam_plot;
figure(1)
subplot(1, 2, 1)
plot(lam_plot / 1e-9, real(er_func(ft_plot * 2 * pi)), 'r-')
xlabel('Wavelength [nm]')
axis tight
subplot(1, 2, 2)
plot(lam_plot / 1e-9, -imag(er_func(ft_plot * 2 * pi)), 'b-')
xlabel('Wavelength [nm]')
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
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 2);
fprintf(fileID, "%d\n", PML_d);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 1 {\n");
% medium 0, scatterer medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 3\n", 1.1431, 0, 1, 0);
fprintf(fileID, "{ Drude %e %e }\n", 1.3202e16/(2*pi), 1.0805e14/(2*pi));
fprintf(fileID, "{ CP %e %e %e %e}\n", 0.26698, -1.2371, 3.8711e15/(2*pi), 4.4642e14/(2*pi));
fprintf(fileID, "{ CP %e %e %e %e}\n", 3.0834, -1.0968, 4.1684e15/(2*pi), 2.3555e15/(2*pi));
fprintf(fileID, "}\n");
% medium end
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
fprintf(fileID, "plane %d %e %e", dim(3), fp, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% scattering coefficients calculated from poynting flux
fprintf(fileID, "flux %s %s\n", flux_input_file, flux_output_file);

disp('config.in created');

% theoretical fields
c_scat_phy = zeros(size(m));

for i = 1 : size(m, 1)
    res = Mie(m(i), size_param(i));
    c_scat_phy(i) = res(2);
end


%% numerical fields
c_scat = load(flux_output_file);
c_scat = c_scat ./ (0.5 * abs(ref).^2 / eta0 * pi * a^2);
disp('Numerical Fields Extracted');

%% comparison
figure(3)
plot(lam / 1e-9, c_scat, 'b*'), hold on
plot(lam / 1e-9, c_scat_phy, 'r-'), hold off
xlabel('Wavelength [nm]')
