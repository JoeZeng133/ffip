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
dim = [52, 52, 50];
time_step = 10000; %gold
PMl_d = 8;

ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);

ft = c0 / (800e-9); %frequency to extract fields at
fp = c0 / (800e-9); %center frequency of the pulse
delay = 0;
t = (0:time_step) * dt;
% ref_signal = ricker(t, fp, delay);
ref_signal = sin(2 * pi * ft * (t - delay));
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * ft * t));
fprintf('Equal to %.2f times of full wavelength\n', time_step * dt * fp);

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
inhom_center = [center(1), center(2), 45e-9 ];
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
sim_x = (-8:60) * dx;
sim_y = (-8:60) * dx;
sim_z = 55e-9;
[sim_X, sim_Y, sim_Z] = ndgrid(sim_x, sim_y, sim_z);
probes_pos = [sim_X(:), sim_Y(:), sim_Z(:)];
num_probes = size(probes_pos, 1);

% convergence check probes
conv_x = (0:4:52) * dx;
conv_y = (0:4:52) * dx;
conv_z = 55e-9;
[conv_X, conv_Y, conv_Z] = ndgrid(conv_x, conv_y, conv_z);
conv_pos = [conv_X(:), conv_Y(:), conv_Z(:)];
num_conv = size(conv_pos, 1);

% write to geometry file
file_geometry = 'vis_geometry.in';
fileID = fopen(file_geometry, 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
for i = 1 : 3
    fprintf(fileID, '%e %e\n', inhom_p1(i), inhom_dx(i));
end
fprintf(fileID, '%e ', rho_target(:));
fclose(fileID);

% write to probes file
file_probes_input = 'vis_probes.in';
fileID = fopen(file_probes_input, "w");
fprintf(fileID, "%d\n", num_probes);
fprintf(fileID, "%e %e %e %e\n", [probes_pos, ft * ones([num_probes 1])]');
fclose(fileID);

% write to convergence check file
file_conv_input = 'conv.in';
fileID = fopen(file_conv_input, "w");
fprintf(fileID, "%d\n", num_conv);
fprintf(fileID, "%e %e %e\n", conv_pos');
fclose(fileID);

disp('geometry, probes files created');
% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 2);
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
Au2_print(fileID);
% medium 2 FePt
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", real(fept_er), fept_sig, ur_bg, 0);
fprintf(fileID, " }\n");

fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 2 {\n");
% geometry 0 gold aperture
fprintf(fileID, "{ ");
fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry);
fprintf(fileID, " }\n");
% geometry 1 Gold layer
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 1, [dx, dx, inhom_p1(3)], [p2(1) - dx, p2(2) - dx, inhom_p2(3)]);
fprintf(fileID, " }\n");
% % geometry 2 FePt Layer
% fprintf(fileID, "{ ");
% fprintf(fileID, "box %d %e %e %e %e %e %e", 2, fept_p1, fept_p2);
% fprintf(fileID, " }\n");

fprintf(fileID, "}\n");

% plane wave source
fprintf(fileID, "source 1 {\n");
fprintf(fileID, "{ ");
fprintf(fileID, "plane %d %e %e", dim(3), fp, delay);
fprintf(fileID, " }\n");
fprintf(fileID, "}\n");

% probes
file_probes_output = 'vis_output.out';
fprintf(fileID, "nearfield %s %s\n", file_probes_input, file_probes_output);

file_conv_output = 'conv.out';
fprintf(fileID, 'nearfield_time %s %s\n', file_conv_input, file_conv_output);
fclose(fileID);

disp('objective configuration created');
%% simulation
call_exe('std_config');

%% load probes data
data = load(file_probes_output);
make_complex = @(x, y) x + 1j * y;

E_vis = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_vis = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E_vis = E_vis ./ ref_signal_fft;
H_vis = H_vis ./ ref_signal_fft;

%% probes visualization
figure
Ex = abs(E_vis(:, 1));
Ey = abs(E_vis(:, 2));
Ez = abs(E_vis(:, 3));
Eabs = abs(Ex.^2 + Ey.^2 + Ez.^2);

F = Ex;
F = reshape(F, size(sim_X));
surf(sim_x, sim_y, F')
xlabel('x')
ylabel('y')

%% convergence data time domain
fileID = fopen(file_conv_output, 'r');
data = fscanf(fileID, '%e');
fclose(fileID);

data = reshape(data, 6, num_conv, []);
figure
for i = 1 : 30 : num_conv
    tmp = data(1, i, :);
    plot(tmp(:)), hold on
end

%% snapshot data loading
fileID = fopen('snapshot.out', 'r');
snapshot_dim = fscanf(fileID, '%d %d %d', [1 3]);
snapshot = fscanf(fileID, '%e');
fclose(fileID);

%% snapshot data visualization
snapshot = reshape(snapshot, snapshot_dim(2), snapshot_dim(3), []); 
norm = max(abs(snapshot(:)));

figure
for i = 1 : 1: size(snapshot, 3)
    ex = snapshot(:, :, i);
    surf(ex / norm * 20)
    axis equal
    zlim([-20 20]);
    getframe;
end

%% projector data visualization
projector = load('projector.out');

norm = max(abs(projector(:)));

for i = 1 : size(projector, 1)
    plot(projector(i, :) / norm)
%     pause(1e-1)
    ylim([-1 1])
    getframe;
end

