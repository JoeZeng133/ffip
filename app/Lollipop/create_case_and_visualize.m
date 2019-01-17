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
dim = [90, 95, 115];
time_step = 10000; %gold
PMl_d = 8;

ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);

ft = c0 / (500e-9); %frequency to extract fields at
fp = c0 / (500e-9); %center frequency of the pulse
delay = 0;
t = (0:time_step) * dt;
% ref_signal = ricker(t, fp, delay);
ref_signal = sin(2 * pi * ft * (t - delay));
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * ft * t));
fprintf('Equal to %.2f times of full wavelength\n', time_step * dt * fp);

center = dim / 2 * dx;
p2 = dim * dx;

disk_center = [90, 95, 95] * 1e-9;
disk_height = 35e-9;
disk_radius = 95e-9;

box_center = [90, 95, 190] * 1e-9;
box_dim = [35, 20, 40] * 1e-9;
box_p1 = box_center - box_dim * 0.5;
box_p2 = box_center + box_dim * 0.5;

% nearfield probe
probe_x = linspace(0, p2(1), 30);
probe_y = linspace(0, p2(2), 30);
probe_z = 210e-9;
[probe_X, probe_Y, probe_Z] = ndgrid(probe_x, probe_y, probe_z);
probes_pos = [probe_X(:), probe_Y(:), probe_Z(:)];
num_probes = size(probes_pos, 1);

% convergence check probes
conv_x = probe_x(1:3:end);
conv_y = probe_y(1:3:end);
conv_z = probe_z;
[conv_X, conv_Y, conv_Z] = ndgrid(conv_x, conv_y, conv_z);
conv_pos = [conv_X(:), conv_Y(:), conv_Z(:)];
num_conv = size(conv_pos, 1);

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

fileID = fopen('config.in', 'w');

% tf padded layer
fprintf(fileID, 'padded_region %d\n', 2);

% basic configuration
fprintf(fileID, "basic {\n");
fprintf(fileID, "%e %e\n", dt, dx);
fprintf(fileID, "%d %d %d\n", dim);
fprintf(fileID, "%d\n", 10);
fprintf(fileID, "%d\n", PMl_d);
fprintf(fileID, "%d\n", time_step);
fprintf(fileID, "%e %e\n", er_bg, ur_bg);
fprintf(fileID, "}\n");

% medium configuration
fprintf(fileID, "medium 2 {\n");
% medium 0, background medium
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
fprintf(fileID, " }\n");
% medium 1 Au
Au_print(fileID);
fprintf(fileID, "}\n");

% geometry configuration
fprintf(fileID, "geometry 2 {\n");
% geometry 0 gold aperture
fprintf(fileID, "{ ");
fprintf(fileID, "disk %d %e %e %e %e %e %d", 1, disk_center, disk_radius, disk_height, 0);
fprintf(fileID, " }\n");
% geometry 1 Gold layer
fprintf(fileID, "{ ");
fprintf(fileID, "box %d %e %e %e %e %e %e", 1, box_p1, box_p2);
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
E = sqrt(Ex.^2 + Ey.^2 + Ez.^2);

F = E;
F = reshape(F, size(probe_X));
pcolor(probe_x/1e-9, probe_y/1e-9, F')
colorbar
shading flat
colormap jet
xlabel('x [nm]')
ylabel('y [nm]')

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

