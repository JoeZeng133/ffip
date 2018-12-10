clear
clc
close all

load('case_configuration.mat');
%% write geometry, probe configurations
rho = ones(size(ch_func));                          %start with rho=1

% write to geometry file
filename_geometry_forward = 'forward_geometry.in';
fileID = fopen(['../', filename_geometry_forward], 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
fprintf(fileID, '%e %e\n', inhom_p1(1) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(2) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(3) * dx, dx);
fprintf(fileID, '%e ', rho(:));
fclose(fileID);

% write to probes file
filename_probes_forward = 'forward_probes.in';
fileID = fopen(['../', filename_probes_forward], "w");
fprintf(fileID, "%d\n", num_probes + num_inhom_pos);
fprintf(fileID, "%e %e %e %e\n", [[probes_pos; inhom_pos] * dx, fp * ones([num_probes + num_inhom_pos, 1])]');
fclose(fileID);

disp('forward simulation geometry and probes files created');
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
fprintf(fileID, "inhom %d %d %s", 1, 0, filename_geometry_forward);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "eigen %d %e %e", dim(3), fp, 0);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
fprintf(fileID, "probe %s\n", filename_probes_forward);
fclose(fileID);

disp('forward simulation config.in created');

%% simulated fields
data = load('../output.out');
make_complex = @(x, y) x + 1j * y;

E_forward = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_forward = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_forward = E_forward / ref_signal_fft;
H_forward = H_forward / ref_signal_fft;

E_probe_forward = E_forward(1:num_probes, :);
H_probe_forward = H_forward(1:num_probes, :);
E_inhom_forward = E_forward(num_probes + 1:end, :);
H_inhom_forward = H_forward(num_probes + 1:end, :);

save("forward_results", "E_probe_forward", "H_probe_forward", "E_inhom_forward", "H_inhom_forward", "rho", "filename_geometry_forward", "filename_probes_forward");

disp('forward simulation fields saved');