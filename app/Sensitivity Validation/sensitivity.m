%% start of optimizaiton
clear
clc
close all

load('target.mat');
er_const = 1;

%%
Omega = 2 * pi * fp;
sig_const = -imag(er_const) * (2 * pi * fp * e0);
er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));

% probes include the original probe and all discretized points in the box
all_pos = [25 25 10; inhom_pos];
file_probes_input_forward = 'forward_probes.in';
fileID = fopen(file_probes_input_forward, "w");
fprintf(fileID, "%d\n", size(all_pos, 1));
fprintf(fileID, "%e %e %e %e\n", [all_pos * dx, fp * ones([size(all_pos, 1), 1])]');
fclose(fileID);

% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 1);
fprintf(fileID, "%d\n", PMl_d);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 1 {\n");
% medium 0, scatterer medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", real(er_const), sig_const, ur_bg, 0);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 1 {\n");

% geometry 0, a box specified by two corner points
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 0, inhom_p1 * dx, inhom_p2 * dx);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "plane %d %e %e", dim(3), fp, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% stop step output
fprintf(fileID, "Stop_Step_Output\n");

% probes
file_output_probes_forward = 'forward_output.out';
fprintf(fileID, "nearfield %s %s\n", file_probes_input_forward, file_output_probes_forward);
fclose(fileID);

disp('forward simulation configuration complete');

%% forward simulation results
call_exe('std_config')
data = load(file_output_probes_forward);
make_complex = @(x, y) x + 1j * y;

E_forward = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_forward = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_forward = E_forward / ref_signal_fft;
H_forward = H_forward / ref_signal_fft;

E_forward_probes = E_forward(1, :);
H_forward_probes = H_forward(1, :);
E_forward_inhom = E_forward(2:end, :);
H_forward_inhom = H_forward(2:end, :);

target_func = 0.5 * sum(abs(E_forward_probes - E_target).^2);
fprintf('Current Target value is %e\n', target_func);
save("forward_results");
disp('forward simulation fields saved');

%% adjoint simulation
Ex = 1;
Ey = 2;
Ez = 4;

%     dipoles = (abs(E_forward_probes) - abs(E_target)) .* conj(E_forward_probes) ./ abs(E_forward_probes);
dipoles = conj(E_forward_probes - E_target);
amp = abs(dipoles);
dipoles_delay = delay - angle(dipoles) / (2 * pi * fp);
output_dipoles = zeros([3 7]);
ctype = [Ex, Ey, Ez];

for i = 1 : 3
    output_dipoles(i, :) = [[25 25 10] * dx, amp(i), fp, dipoles_delay(i), ctype(i)];
end

file_dipoles_adjoint = 'adjoint_dipoles.in';
fileID = fopen(file_dipoles_adjoint, 'w');
fprintf(fileID, '%d\n', size(output_dipoles, 1));
fprintf(fileID, '%e %e %e %e %e %e %d\n', output_dipoles');
fclose(fileID);

% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 1);
fprintf(fileID, "%d\n", PMl_d);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 1 {\n");
% medium 0, scatterer medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", real(er_const), sig_const, ur_bg, 0);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 1 {\n");
% geometry 0, a box specified by two corner points
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 0, inhom_p1 * dx, inhom_p2 * dx);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% dipole sources
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "dipole %s", file_dipoles_adjoint);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% stop step output
fprintf(fileID, "Stop_Step_Output\n");

% probes
file_probes_output_adjoint = 'adjoint_output.out';
fprintf(fileID, "nearfield %s %s\n", file_probes_input_forward, file_probes_output_adjoint );
fclose(fileID);

disp('adjoint simulation configuration complete');
%% adjoint fields
call_exe('std_config');
data = load(file_probes_output_adjoint);
E_adjoint = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_adjoint = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_adjoint = E_adjoint / ref_signal_fft;
H_adjoint = H_adjoint / ref_signal_fft;

E_adjoint_probes = E_adjoint(1, :);
H_adjoint_probes = H_adjoint(1, :);
E_adjoint_inhom = E_adjoint(2:end, :);
H_adjoint_inhom = H_adjoint(2:end, :);

save("adjoint_results");
disp('adjoint simulation results saved');
% Sensitivity Calculation
int_weights = @(dim) [0.5, ones([1, dim - 2]), 0.5];
[w_X, w_Y, w_Z] = ndgrid(int_weights(inhom_dim(1) + 1), int_weights(inhom_dim(2) + 1), int_weights(inhom_dim(3) + 1));
w = w_X .* w_Y .* w_Z;
w = w(:) * dx^3;
Se = 1j * (2 * pi * fp) * sum(E_forward_inhom .* E_adjoint_inhom, 2);
Se = sum(Se .* w) * e0;
fprintf('Se = ')
disp(Se);

% Calculate updates
step_er = 0.05;
f = [real(Se); -imag(Se)];
lb = [max(1 - real(er_const), -step_er); max(-1 - imag(er_const), -step_er)];
ub = [min(3 - real(er_const), step_er); min(1 - imag(er_const), step_er)];
x = linprog(f, [], [], [], [], lb, ub);

delta_er = x(1) + 1j * x(2);
new_target_func = real(Se * delta_er) + target_func;
fprintf('new er is going to be :');
disp(er_const + delta_er);
fprintf('Expected Target Function Value is going to be : %e\n', new_target_func);

% Update er
er_const = er_const + delta_er;
disp('er updated');





