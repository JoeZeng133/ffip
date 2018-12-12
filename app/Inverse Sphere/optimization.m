clear
clc
close all

load('case_configuration.mat');
%% forward simulation
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

%% forward simulation results
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

%% adjoint simulation
Ex = 1;
Ey = 2;
Ez = 4;
dipoles = 2 * (abs(E_forward_probes) - abs(E_target_probes)) .* (E_forward_probes) ./ abs(E_forward_probes);
amp = abs(dipoles);
delay = (angle(dipoles) + pi) / (2 * pi * fp) + d;
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

disp('adjoint config.in created');
%% adjoint fields
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
%% Sensitivity Calculation
int_weights = @(dim) [0.5, ones([1, dim - 2]), 0.5];
[w_X, w_Y, w_Z] = ndgrid(int_weights(inhom_dim(1) + 1), int_weights(inhom_dim(2) + 1), int_weights(inhom_dim(3) + 1));
w = w_X .* w_Y .* w_Z;
w = w(:) * dx^3;
Se = 1j * (2 * pi * fp) * sum(E_forward_inhom .* E_adjoint_inhom, 2);
drho = e0 * (er_const - er_bg);
A = real(Se .* w * drho);
A = reshape(A, inhom_dim + 1);
target_func = sum(sum((abs(E_forward_probes) - abs(E_target_probes)).^2, 2), 1);

%% Visualizing Sensitivity
for i = 1 : 4
    subplot(2, 2, i)
    level = floor(i * (inhom_dim(3) - 1) / 4 + 1);
    surf(A(:, :, level));
    xlabel('y')
    ylabel('x');
    title(num2str(level))
end

