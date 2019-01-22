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
dim = [50, 50, 50];
time_step = 10000;
tf_layer = 2;
sf_layer = 1;
projector_padding = ceil((tf_layer + 1) / 2);
PML_d = 8;

ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
fs = c0 / (800e-9); %center frequency of the pulse
ft = c0 / (800e-9); %frequency to extract fields at
delay = 0;
fprintf('Equal to %.2f times of full pulse width\n', time_step * dt * fs);

center = dim / 2 * dx;
p2 = dim * dx;

% FePt 
fept_er = (2.9 - 1.5j)^2;
fept_sig = -imag(fept_er) * (2 * pi * ft * e0);
fept_len = 100e-9;
fept_h = 10e-9;
fept_center = [center(1), center(2), 60e-9];
fept_p1 = fept_center - [fept_len, fept_len, fept_h] / 2;
fept_p2 = fept_center + [fept_len, fept_len, fept_h] / 2;

er_const = Au(2 * pi * ft); %gold

% disk aperture
disk_length = 70e-9;
disk_radius = 30e-9;
disk_h = 10e-9;
inhom_center = [center(1), center(2), 45e-9];
inhom_p1 = inhom_center - [disk_length, disk_length, disk_h] / 2;
inhom_p2 = inhom_center + [disk_length, disk_length, disk_h] / 2;
inhom_dim = ceil((inhom_p2 - inhom_p1) / dx);
inhom_dx = (inhom_p2 - inhom_p1) ./ inhom_dim;
inhom_coord = cell([1 3]);
inhom_coord_grid = cell([1 3]);
for i = 1 : 3
    inhom_coord{i} = linspace(inhom_p1(i), inhom_p2(i), inhom_dim(i) + 1);
end
[inhom_coord_grid{1}, inhom_coord_grid{2}, inhom_coord_grid{3}] = ndgrid(inhom_coord{1}, inhom_coord{2}, inhom_coord{3});
inhom_pos = [inhom_coord_grid{1}(:), inhom_coord_grid{2}(:), inhom_coord_grid{3}(:)];
num_inhom_pos = size(inhom_pos, 1);
rho_target = ((inhom_coord_grid{1} - center(1)).^2 + (inhom_coord_grid{2} - center(2)).^2 > disk_radius^2) * 1;

% nearfield probe
target_x = linspace(0, p2(1), 10);
target_y = linspace(0, p2(2), 10);
target_z = 60e-9;
[target_X, target_Y, target_Z] = ndgrid(target_x, target_y, target_z);
probes_pos = [target_X(:), target_Y(:), target_Z(:)];
num_probes = size(probes_pos, 1);

% write to geometry file
file_geometry_target = 'geometry_target.in';
fileID = fopen(file_geometry_target, 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
for i = 1 : 3
    fprintf(fileID, '%e %e\n', inhom_p1(i), inhom_dx(i));
end
fprintf(fileID, '%e\n', rho_target(:));
fclose(fileID);

% write to probes file
file_probes_input_target = 'probes.in';
fileID = fopen(file_probes_input_target, "w");
fprintf(fileID, "%d\n", num_probes);
fprintf(fileID, "%e %e %e %e\n", [probes_pos, ft * ones([num_probes 1])]');
fclose(fileID);

disp('geometry, probes files created');
% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, 'basic {\n');
fprintf(fileID, '%e %e\n', dt, dx);
fprintf(fileID, '%d %d %d\n', dim);
fprintf(fileID, '%d\n', sf_layer);
fprintf(fileID, '%d\n', tf_layer);
fprintf(fileID, '%d\n', PML_d);
fprintf(fileID, '%d\n', time_step);
fprintf(fileID, '%e %e\n', er_bg, ur_bg);
fprintf(fileID, '}\n');

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
fprintf(fileID, "%e %e %e %e 0", real(fept_er), fept_sig, ur_bg, 0);
fprintf(fileID, " }\n");

fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 3 {\n");
% geometry 0 gold aperture
fprintf(fileID, "{ ");
fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_target);
fprintf(fileID, " }\n");
% geometry 1 Gold layer
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 1, [0, 0, inhom_p1(3)], [p2(1), p2(2), inhom_p2(3)]);
fprintf(fileID, " }\n");
% % geometry 2 FePt Layer
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 2, fept_p1, fept_p2);
fprintf(fileID, " }\n");

fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, 'plane %d %d %s %e %e %e', projector_padding, dim(3) + projector_padding, 's', fs, delay, 0);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
file_probes_output_target = 'target_output.out';
fprintf(fileID, "nearfield %s %s\n", file_probes_input_target, file_probes_output_target);

disp('objective configuration created');
%% simulated fields
call_exe('std_config');

%% load data
ref_signal = load('reference.out');
ref_signal = reshape(ref_signal, 1, []);
t = (0:numel(ref_signal) - 1) * dt;
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * ft(:) * t), 2);

data = load(file_probes_output_target);
make_complex = @(x, y) x + 1j * y;

E_target_probes = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_target_probes = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E_target_probes = E_target_probes / ref_signal_fft;
H_target_probes = H_target_probes / ref_signal_fft;

save('case_configuration');
disp('objective fields saved');

%% visualization
figure
v_rho_target = reshape(rho_target, inhom_dim + 1);
v_rho_target = v_rho_target(:, :, 1);
pcolor(v_rho_target)
colorbar
shading flat
title('The original')
xlabel('y')
ylabel('x')

figure
H = sqrt(sum(abs(H_target_probes).^2, 2));
H = reshape(H, size(target_X));
surf(target_x, target_y, H')
shading flat
title('|H|')

figure
E = sqrt(sum(abs(E_target_probes).^2, 2));
E = reshape(E, size(target_X));
surf(target_x, target_y, E')
shading flat
title('|E|')
