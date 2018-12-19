clear
clc
close all

%% case configuration
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
Np = 30;                           
fp = c0 / (Np * dx);
delay = 1 / fp;
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, fp, delay);
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * fp * t));

% medium for the box
Omega = 2 * pi * fp;
er_const = 1.2 - 0.2i;
sig_const = -imag(er_const) * (2 * pi * fp * e0);
er_func = @(w) (real(er_const) - 1j * sig_const ./ (w * e0));

% Box geometry [15, 15, 15] to [35, 35, 35] in [50, 50, 50]
inhom_p1 = [15, 15, 15];
inhom_dim = [20, 20, 20];
inhom_p2 = inhom_dim + inhom_p1;
[inhom_x, inhom_y, inhom_z] = ndgrid(inhom_p1(1):inhom_p2(1), inhom_p1(2):inhom_p2(2), inhom_p1(3):inhom_p2(3));
inhom_pos = [inhom_x(:), inhom_y(:), inhom_z(:)];

% put one probe at [25, 25, 10]
all_pos = [25 25 10];
file_probes_input_target = 'target_probes.in';
fileID = fopen(file_probes_input_target, "w");
fprintf(fileID, "%d\n", 1);
fprintf(fileID, "%e %e %e %e\n", [all_pos * dx, fp]);
fclose(fileID);

% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "%d\n", PMl_d);
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
fprintf(fileID, "eigen %d %e %e", dim(3), fp, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
file_output_probes_target = 'target_output.out';
fprintf(fileID, "probe %s %s\n", file_probes_input_target, file_output_probes_target);
fclose(fileID);

disp('case configuration done');
%% forward simulation results
!./std_config
data = load(file_output_probes_target);
make_complex = @(x, y) x + 1j * y;

E_target = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_target = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_target = E_target / ref_signal_fft;
H_target = H_target / ref_signal_fft;

disp('case simulation done');
save('target');