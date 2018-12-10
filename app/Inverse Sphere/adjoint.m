clear
clc
close all

load('case_configuration.mat');
load('forward_results.mat');

%% write dipoles
Ex = 1;
Ey = 2;
Ez = 4;
dipoles = 2 * (abs(E_probe_forward) - abs(E_probe_objective)) .* conj(E_probe_forward) ./ abs(E_probe_forward);
amp = abs(dipoles);
delay = angle(dipoles) / (2 * pi * fp);
output_dipoles = [probes_pos * dx, amp(:, 1), fp * ones([num_probes, 1]), delay(:, 1), Ex * ones([num_probes, 1]);...
    probes_pos * dx, amp(:, 1), fp * ones([num_probes, 1]), delay(:, 1), Ey * ones([num_probes, 1]);...
    probes_pos * dx, amp(:, 1), fp * ones([num_probes, 1]), delay(:, 1), Ez * ones([num_probes, 1]);];

filename_dipoles = 'dipoles.in';
fileID = fopen(['../', filename_dipoles], 'w');
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
fprintf(fileID, "inhom %d %d %s", 1, 0, filename_geometry_forward);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% dipole sources
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "dipole %s", filename_dipoles);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
fprintf(fileID, "probe %s\n", filename_probes_forward);
fclose(fileID);

disp('adjoint simulation config.in creatd');

%% simulated fields
data = load('../output.out');
make_complex = @(x, y) x + 1j * y;

E_adjoint = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_adjoint = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_adjoint = E_adjoint / ref_signal_fft;
H_adjoint = H_adjoint / ref_signal_fft;

E_probe_adjoint = E_adjoint(1:num_probes, :);
H_probe_adjoint = H_adjoint(1:num_probes, :);
E_inhom_adjoint = E_adjoint(num_probes + 1:end, :);
H_inhom_adjoint = H_adjoint(num_probes + 1:end, :);

save("adjoint_results", "E_probe_adjoint", "H_probe_adjoint", "E_inhom_adjoint", "H_inhom_adjoint");
