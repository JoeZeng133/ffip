clear
clc
close all

load('forward_results.mat');

%% write dipoles
Ex = 1;
Ey = 2;
Ez = 4;
dipoles = 2 * (abs(E_forward_probes) - abs(E_target_probes)) .* conj(E_forward_probes) ./ abs(E_forward_probes);
amp = abs(dipoles);
delay = angle(dipoles) / (2 * pi * fp) + d;
output_dipoles = [...
    probes_pos * dx, amp(:, 1), fp * ones([num_probes, 1]), delay(:, 1), Ex * ones([num_probes, 1]);...
    probes_pos * dx, amp(:, 2), fp * ones([num_probes, 1]), delay(:, 2), Ey * ones([num_probes, 1]);...
    probes_pos * dx, amp(:, 3), fp * ones([num_probes, 1]), delay(:, 3), Ez * ones([num_probes, 1]);];

file_dipoles_adjoint = 'adjoint_dipoles.in';
fileID = fopen(['../', file_dipoles_adjoint], 'w');
fprintf(fileID, '%d\n', num_probes * 3);
fprintf(fileID, '%e %e %e %e %e %e %d\n', output_dipoles');
fclose(fileID);

disp('adjoint simulation dipoles file created');

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

% dipole sources
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "dipole %s", file_dipoles_adjoint);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
file_probes_output_adjoint = 'adjoint_output.out';
fprintf(fileID, "probe %s %s\n", file_probes_input_forward, file_probes_output_adjoint );
fclose(fileID);

disp('config.in created');

%% simulated fields
data = load(['../', file_probes_output_adjoint ]);

E_adjoint = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_adjoint = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_adjoint = E_adjoint / ref_signal_fft;
H_adjoint = H_adjoint / ref_signal_fft;

E_adjoint_probes = E_adjoint(1:num_probes, :);
H_adjoint_probes = H_adjoint(1:num_probes, :);
E_adjoint_inhom = E_adjoint(num_probes + 1:end, :);
H_adjoint_inhom = H_adjoint(num_probes + 1:end, :);

save("adjoint_results");
