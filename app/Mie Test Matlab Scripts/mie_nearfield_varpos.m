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
sf_layer = 1;
projector_padding = ceil((tf_layer + 1) / 2);

fs = c0 / 500e-9;
delay = 1 / fs;

lambda = 800e-9;
ft = c0 / lambda;
a = 40e-9;

center = dim * dx / 2;
rho = a + linspace(10e-9, 20e-9, 10);
phi = linspace(0, 2 * pi, 10);
th = linspace(0,  pi, 10);
[xc, yc, zc] = sph2cart(phi, pi / 2 - th, rho);

[Xc, Yc, Zc, Ft] = ndgrid(xc, yc, zc, ft);
probes = [Xc(:) + center(1), Yc(:) + center(2), Zc(:) + center(3), Ft(:)];

% near field probes
filename_probes = 'probes.in';
fileID = fopen(filename_probes, 'w');
fprintf(fileID, '%d\n', size(probes, 1));
fprintf(fileID, '%e %e %e %e\n', probes');
fclose(fileID);

Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @Au;
er = er_func(Omega);

ns = sqrt(conj(er_func(2 * pi * ft)));
nm = 1;

% figure(1)
% ft_samples = linspace(c0 / 1000e-9, c0 / 200e-9, 100);
% plot(c0 ./ ft_samples / 1e-9, real(er_func(ft_samples * 2 * pi)), 'r-'), hold on
% plot(c0 ./ ft_samples / 1e-9, -imag(er_func(ft_samples * 2 * pi)), 'b-')
% xlabel('Wavelength [nm]')
% legend({'$\textrm{Re}(\varepsilon_r)$', '$\textrm{Im}(\varepsilon_r)$'}, 'interpreter', 'latex','fontsize', 15)
% axis tight

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
fprintf(fileID, 'plane %d %d %s %e %e %e', projector_padding, dim(3) + projector_padding, 'r', fs, delay, dim(3) / 2 * dx);
fprintf(fileID, ' }\n');
fprintf(fileID, '}\n');

% farfield probes
fprintf(fileID, 'nearfield %s %s\n', filename_probes, 'output.out');
fclose(fileID);

disp('config.in created');

%% theoretical fields
[E, H] = calcmie_nf( a, ns, nm, lambda, Xc(:), Yc(:), Zc(:), 'TotalField', true );

%% numerical fields
data = load('output.out');
make_complex = @(x, y) x + 1j * y;
ref_signal = load('reference.out');
ref_signal = reshape(ref_signal, 1, []);
t = (0:numel(ref_signal) - 1) * dt;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

En = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
Hn = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
En = En ./ ref;
Hn = Hn ./ ref;

disp('Numerical Fields Extracted');


%% corellation plots
for i = 1 : 3
    figure(i + 1)
    subplot(2, 2, 1)
    plot(real(E(:, i)), real(En(:, i)), '.'), hold on
    plot(real(E(:, i)), real(E(:, i))), hold off
    xlabel('Analytical')
    axis equal
    
    subplot(2, 2, 2)
    plot(imag(E(:, i)), -imag(En(:, i)), '.'), hold on
    plot(imag(E(:, i)), imag(E(:, i))), hold off
    xlabel('Analytical')
    axis equal

    subplot(2, 2, 3)
    plot(real(H(:, i)), real(Hn(:, i)), '.'), hold on
    plot(real(H(:, i)), real(H(:, i))), hold off
    xlabel('Analytical')
    axis equal
    
    subplot(2, 2, 4)
    plot(imag(H(:, i)), -imag(Hn(:, i)), '.'), hold on
    plot(imag(H(:, i)), imag(H(:, i))), hold off
    xlabel('Analytical')
    axis equal
end


