clear
clc
close all

%% configuration input files
% constants
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

% basic configuration
er_bg = 1;
ur_bg = 1;
Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [50, 50, 50];
time_step = 600;
PMl_d = 6;

% excitation reference and frequency of interests (1 frequency)
Np = 30;                           
fp = c0 / (Np * dx);
delay = 1 / fp;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:time_step) * dt;
ref_signal = ricker(t, fp, delay);
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * fp * t));

% medium for the inhomogeneous region
Omega = 2 * pi * fp;
er_const = 2 - 0.2i;
sig_const = -imag(er_const) * (2 * pi * fp * e0);
er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));

% Geometry for the inhomogeneous region
sphere_func = @(X, Y, Z, r, center) ((X - center(1)).^2 + (Y - center(2)).^2) <= r^2;
geo_func = @(X, Y, Z) sphere_func(X, Y, Z, 10, dim / 2);
inhom_p1 = [10, 10, 20];
inhom_dim = [30, 30, 10];
inhom_p2 = inhom_p1 + inhom_dim;
inhom_x = inhom_p1(1):inhom_p2(1);
inhom_y = inhom_p1(2):inhom_p2(2);
inhom_z = inhom_p1(3):inhom_p2(3);
[inhom_X, inhom_Y, inhom_Z] = ndgrid(inhom_x, inhom_y, inhom_z);
inhom_pos = [inhom_X(:), inhom_Y(:), inhom_Z(:)];
num_inhom_pos = size(inhom_pos, 1);
rho_target = geo_func(inhom_X, inhom_Y, inhom_Z) * 1;

% num_probes = 10;
% [sim_X, sim_Y, sim_Z] = ndgrid(0:dim(1), 0:dim(2), 0:dim(3));
% slice = ...
%     (sim_X < inhom_p1(1) - 5 | sim_X > inhom_p2(1)) + 5 | ...
%     (sim_Y < inhom_p1(2) - 5 | sim_Y > inhom_p2(2)) + 5 | ...
%     (sim_Z < inhom_p1(3) - 5 | sim_Z > inhom_p2(3)) + 5;
% 
% sim_X = sim_X(slice);
% sim_Y = sim_Y(slice);
% sim_Z = sim_Z(slice);
% probes_pos = datasample([sim_X(:), sim_Y(:), sim_Z(:)], num_probes, 1, 'Replace', false);

[sim_X, sim_Y, sim_Z] = ndgrid([10, 20, 30, 40], [10, 20, 30, 40], [5, 10]);
probes_pos = [sim_X(:), sim_Y(:), sim_Z(:)];
num_probes = size(probes_pos, 1);

% write to geometry file
file_geometry_target = 'target_geometry.in';
fileID = fopen(file_geometry_target, 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
fprintf(fileID, '%e %e\n', inhom_p1(1) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(2) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(3) * dx, dx);
fprintf(fileID, '%e ', rho_target(:));
fclose(fileID);

% write to probes file
file_probes_input_target = 'target_probes.in';
fileID = fopen(file_probes_input_target, "w");
fprintf(fileID, "%d\n", num_probes);
fprintf(fileID, "%e %e %e %e\n", [probes_pos * dx, fp * ones([num_probes 1])]');
fclose(fileID);

disp('geometry, probes files created');
% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", time_step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "%d\n", PMl_d);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 2 {\n");
% medium 0, background medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
fprintf(fileID, " }\n");

% medium 1, scatterer medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", real(er_const), sig_const, ur_bg, 0);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% geometry configuration

fprintf(fileID, "geometry 1 {\n");
% geometry 0, the inhomogeneous region with mixed medium1 and medium0
fprintf(fileID, "{ ");
fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_target);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "eigen %d %e %e", dim(3), fp, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
file_probes_output_target = 'target_output.out';
fprintf(fileID, "probe %s %s\n", file_probes_input_target, file_probes_output_target);
fclose(fileID);

disp('objective configuration created');
%% simulated fields
!std_config.exe
data = load(file_probes_output_target);
make_complex = @(x, y) x + 1j * y;

E_target_probes = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_target_probes = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E_target_probes = E_target_probes / ref_signal_fft;
H_target_probes = H_target_probes / ref_signal_fft;

save('case_configuration');
disp('objective fields saved');

figure
v_rho_target = reshape(rho_target, inhom_dim + 1);
v_rho_target = v_rho_target(:, :, 1);
pcolor(v_rho_target)
colorbar
shading flat
title('The original')
xlabel('y')
ylabel('x')
