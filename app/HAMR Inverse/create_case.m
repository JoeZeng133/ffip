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
Sc = 0.5;
dx = 2e-9;
dt = Sc * dx / c0;
dim = [50, 50, 30];
time_step = 20000;
PMl_d = 6;

% excitation reference and frequency of interests (1 frequency)                         
fp = c0 / (800e-9);
delay = 1 / fp;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:time_step) * dt;
ref_signal = ricker(t, fp, delay);
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * fp * t));

% medium for FePT layer
Omega = 2 * pi * fp;
fept_er = (2.9 - 1.5j)^2;
fept_sig = -imag(fept_er) * (2 * pi * fp * e0);
er_func = @(w) (real(fept_er) - 1j * fept_sig ./ (w * e0));
fept_p1 = [5e-9, 5e-9, 34e-9];
fept_p2 = [95e-9, 95e-9, 54e-9];

% Geometry for the gold layer
sphere_func = @(X, Y, Z, r, center) ((X - center(1)).^2 + (Y - center(2)).^2) <= r^2;
aperture_p1 = [5e-9, 5e-9, 4e-9];
aperture_p2 = [95e-9,95e-9,24e-9];
aperture_dim = ceil((aperture_p2 - aperture_p1) / dx);
aperture_dx = (aperture_p2 - aperture_p1) ./ aperture_dim;
aperture_points = cell([1 3]);
for i = 1 : 3
    aperture_points{i} = linspace(aperture_p1(i), aperture_p2(i), aperture_dim(i) + 1);
end
[aperture_X, aperture_Y, aperture_Z] = ndgrid(aperture_points{1}, aperture_points{2}, aperture_points{3});
aperture_pos = [aperture_X(:), aperture_Y(:), aperture_Z(:)];
num_aperture_pos = size(aperture_pos, 1);

geo_func = @(X, Y, Z) 1 - sphere_func(X, Y, Z, 30e-9, (aperture_p1 + aperture_p2) / 2);
rho_target = geo_func(aperture_X, aperture_Y, aperture_Z);

% nearfield probe positions
probe_x = linspace(20e-9, 80e-9, 5);
probe_y = linspace(20e-9, 80e-9, 5);
probe_z = 44e-9;
[probe_X, probe_Y, probe_Z] = ndgrid(probe_x, probe_y, probe_z);
probes_pos = [probe_X(:), probe_Y(:), probe_Z(:)];
num_probes = size(probes_pos, 1);

% write to geometry file
file_geometry_target = 'target_geometry.in';
fileID = fopen(file_geometry_target, 'w');
fprintf(fileID, '%d %d %d\n', aperture_dim + 1);
fprintf(fileID, '%e %e\n', aperture_p1(1), aperture_dx(1));
fprintf(fileID, '%e %e\n', aperture_p1(2), aperture_dx(2));
fprintf(fileID, '%e %e\n', aperture_p1(3), aperture_dx(3));
fprintf(fileID, '%e ', rho_target(:));
fclose(fileID);

% write to probes file
file_probes_input_target = 'target_probes.in';
fileID = fopen(file_probes_input_target, "w");
fprintf(fileID, "%d\n", num_probes);
fprintf(fileID, "%e %e %e %e\n", [probes_pos, fp * ones([num_probes 1])]');
fclose(fileID);

disp('geometry, probes files created');
% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 1);
fprintf(fileID, "%d\n", PMl_d);
fprintf(fileID, "%d\n", time_step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 3 {\n");
% medium 0, background medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
fprintf(fileID, " }\n");
% medium 1 Au
Au_print(fileID);
% medium 2 FePt
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0\n", real(fept_er), fept_sig, ur_bg, 0);
fprintf(fileID, "}\n");

fprintf(fileID, " }\n");

% geometry configuration
fprintf(fileID, "geometry 2 {\n");
% geometry 0, the aperture
fprintf(fileID, "{ ");
fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_target);
fprintf(fileID, " }\n");
% geometry 1, the FePt layer 20nm thick
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 2, fept_p1, fept_p2);
fprintf(fileID, " }\n");

fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "plane %d %e %e", dim(3), fp, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% no step number output
% fprintf(fileID, 'Stop_Step_Output\n');

% probes
file_probes_output_target = 'target_output.out';
fprintf(fileID, "nearfield %s %s\n", file_probes_input_target, file_probes_output_target);
fclose(fileID);

disp('objective configuration created');
%% simulated fields
% call_exe('std_config');
data = load(file_probes_output_target);
make_complex = @(x, y) x + 1j * y;

E_target_probes = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_target_probes = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E_target_probes = E_target_probes / ref_signal_fft;
H_target_probes = H_target_probes / ref_signal_fft;

save('case_configuration');
disp('objective fields saved');

figure
v_rho_target = reshape(rho_target, aperture_dim + 1);
v_rho_target = v_rho_target(:, :, 1);
pcolor(v_rho_target)
colorbar
shading flat
title('The original')
xlabel('y')
ylabel('x')
