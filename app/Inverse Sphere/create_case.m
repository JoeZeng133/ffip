clear
clc
close all

%%
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
step = 600;
PMl_d = 6;

% excitation reference and frequency of interests (1 frequency)
Np = 20;                           
fp = c0 / (Np * dx);
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, fp, 0);
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * fp * t), 2);

% medium for the inhomogeneous region
Omega = 2 * pi * fp;
lorentz = @(w, rel_e, fp, delta) rel_e ./ (1 + 1j * (w / (2 * pi * fp)) * (delta / fp) - (w / (2 * pi * fp)).^2);
drude = @(w, fp, gamma) (2 * pi * fp)^2 ./ (1j * w * (2 * pi * gamma) - w.^2);
deybe = @(w, rel_e, tau) rel_e ./ (1 + 1j * w * tau);

er_const = 2 - 0.2i;
sig_const = -imag(er_const) * (2 * pi * fp * e0);
er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));

% Geometry for the inhomogeneous region
sphere_func = @(X, Y, Z, r, center) ((X - center(1)).^2 + (Y - center(2)).^2 + (Z - center(3)).^2) <= r^2;
geo_func = @(X, Y, Z) sphere_func(X, Y, Z, 10, dim / 2);

inhom_p1 = [15, 15, 15];
inhom_dim = [20, 20, 20];
inhom_x = inhom_p1(1) + (0:inhom_dim(1));
inhom_y = inhom_p1(2) + (0:inhom_dim(2));
inhom_z = inhom_p1(3) + (0:inhom_dim(3));
[inhom_X, inhom_Y, inhom_Z] = ndgrid(inhom_x, inhom_y, inhom_z);
inhom_pos = [inhom_X(:), inhom_Y(:), inhom_Z(:)];
num_inhom_pos = size(inhom_pos, 1);
ch_func = geo_func(inhom_X, inhom_Y, inhom_Z);

disp('basic configuring done');
%% probes, outside of inhom region
num_probes = 100;
[sim_X, sim_Y, sim_Z] = ndgrid(0:dim(1), 0:dim(2), 0:dim(3));
slice = ~geo_func(sim_X, sim_Y, sim_Z);
sim_X = sim_X(slice);
sim_Y = sim_Y(slice);
sim_Z = sim_Z(slice);
probes_pos = datasample([sim_X(:), sim_Y(:), sim_Z(:)], num_probes, 1, 'Replace', false);

disp('probes sampling done');

%% write to geometry files
% write to geometry file
filename_geometry_objective = 'objective_geometry.in';
fileID = fopen(['../', filename_geometry_objective], 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
fprintf(fileID, '%e %e\n', inhom_p1(1) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(2) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(3) * dx, dx);
fprintf(fileID, '%e ', ch_func(:));
fclose(fileID);

% write to probes file
filename_probes = 'probes.in';
fileID = fopen(['../', filename_probes], "w");
fprintf(fileID, "%d\n", num_probes);
fprintf(fileID, "%e %e %e %e\n", [probes_pos * dx, fp * ones([num_probes 1])]');
fclose(fileID);

disp('geometry, probes files created');
%% write configurations
% basic configuration
fileID = fopen('../config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", step);
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
fprintf(fileID, "inhom %d %d %s", 1, 0, filename_geometry_objective);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "eigen %d %e %e", dim(3), fp, 0);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
fprintf(fileID, "probe %s\n", filename_probes);
fclose(fileID);

disp('config.in created');

%% simulated fields
data = load('../output.out');
make_complex = @(x, y) x + 1j * y;

E_probe_objective = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_probe_objective = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E_probe_objective = E_probe_objective / ref_signal_fft;
H_probe_objective = H_probe_objective / ref_signal_fft;

save('case_configuration');
disp('objective fields saved');
