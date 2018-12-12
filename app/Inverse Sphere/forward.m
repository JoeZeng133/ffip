clear
clc
close all

load('case_configuration.mat');
%% write geometry, probe configurations
% write to geometry file
file_geometry_forward = 'forward_geometry.in';
fileID = fopen(['../', file_geometry_forward], 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
fprintf(fileID, '%e %e\n', inhom_p1(1) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(2) * dx, dx);
fprintf(fileID, '%e %e\n', inhom_p1(3) * dx, dx);
fprintf(fileID, '%e ', rho(:));
fclose(fileID);

% write to probes file
file_probes_input_forward = 'forward_probes.in';
fileID = fopen(['../', file_probes_input_forward], "w");
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
fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_forward);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "eigen %d %e %e", dim(3), fp, d);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
file_output_probes_forward = 'forward_output.out';
fprintf(fileID, "probe %s %s\n", file_probes_input_forward, file_output_probes_forward);
fclose(fileID);

disp('config.in created');

%% simulated fields
data = load(['../', file_output_probes_forward]);

E_forward = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_forward = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_forward = E_forward / ref_signal_fft;
H_forward = H_forward / ref_signal_fft;

E_forward_probes = E_forward(1:num_probes, :);
H_forward_probes = H_forward(1:num_probes, :);
E_forward_inhom = E_forward(num_probes + 1:end, :);
H_forward_inhom = H_forward(num_probes + 1:end, :);

save("forward_results");
disp('forward simulation fields saved');