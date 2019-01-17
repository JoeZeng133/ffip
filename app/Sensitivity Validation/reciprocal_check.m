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
ref_signal = sin(2 * pi * ft * t);
ref_signal_fft = sum(ref_signal .* exp(-1j * 2 * pi * ft * t));
fprintf('Equal to %.2f times of full pulse width\n', time_step * dt * fp);

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

er_const = Au2(2 * pi * ft); %gold

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

% two dipole positions
G1 = [2 + 5j, 0.5 - 3j, 3 - 5j];
x1 = [30e-9, 47e-9, 45e-9];
G2 = [1 + 1j, 2 - 1j, 1 - 3j];
x2 = [10e-9, 90e-9, 60e-9];

%% E2(x1) * G1
% write to probes file
file_probes_input = 'probes.in';
fileID = fopen(file_probes_input, 'w');
fprintf(fileID, '1\n');
fprintf(fileID, '%e %e %e %e\n', x1, fp);
fclose(fileID);

% write to dipoles file
ctype = [1, 2, 4];
file_dipoles_input = 'dipoles.in';
fileID = fopen(file_dipoles_input, 'w');
fprintf(fileID, '3\n');
for i = 1 : 3
    fprintf(fileID, '%e %e %e %e %e %e %d\n', x2, abs(G2(i)), fp, -angle(G2(i)) / (2 * pi * fp), ctype(i));
end
fclose(fileID);

%% configurations
% write to geometry file
file_geometry = 'geometry.in';
fileID = fopen(file_geometry, 'w');
fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
for i = 1 : 3
    fprintf(fileID, '%e %e\n', inhom_p1(i), inhom_dx(i));
end
fprintf(fileID, '%e ', rho_target(:));
fclose(fileID);

disp('geometry, probes files created');
% basic configuration
fileID = fopen('config.in', 'w');
fprintf(fileID, 'basic {\n');
fprintf(fileID, '%e %e\n', dt, dx);
fprintf(fileID, '%d %d %d\n', dim);
fprintf(fileID, '%d\n', 0);
fprintf(fileID, '%d\n', PMl_d);
fprintf(fileID, '%d\n', time_step);
fprintf(fileID, '%e %e\n', er_bg, ur_bg);
fprintf(fileID, '}\n');

% medium configuration
fprintf(fileID, 'medium 3 {\n');
% medium 0, background medium
fprintf(fileID, '{ ');
fprintf(fileID, '%e %e %e %e 0', er_bg, 0, ur_bg, 0);
fprintf(fileID, ' }\n');
% medium 1 Au
Au2_print(fileID);
% medium 2 FePt
fprintf(fileID, '{ ');
fprintf(fileID, '%e %e %e %e 0', real(fept_er), fept_sig, ur_bg, 0);
fprintf(fileID, ' }\n');

fprintf(fileID, '}\n');

% geometry configuration
fprintf(fileID, 'geometry 2 {\n');
% geometry 0 gold aperture
% fprintf(fileID, '{ ');
% fprintf(fileID, 'inhom %d %d %s', 1, 0, file_geometry);
% fprintf(fileID, ' }\n');
% geometry 1 Gold layer
fprintf(fileID, '{ ');
fprintf(fileID, 'box %d %e %e %e %e %e %e', 1, [dx, dx, inhom_p1(3)], [p2(1) - dx, p2(2) - dx, inhom_p2(3)]);
fprintf(fileID, ' }\n');
% geometry 2 FePt Layer
fprintf(fileID, '{ ');
fprintf(fileID, 'box %d %e %e %e %e %e %e', 2, fept_p1, fept_p2);
fprintf(fileID, ' }\n');

fprintf(fileID, '}\n');

fprintf(fileID, 'Stop_Step_Output\n');

% plane wave source
fprintf(fileID, 'source 1 {\n');
fprintf(fileID, '{ ');
fprintf(fileID, 'dipole %s', file_dipoles_input);
fprintf(fileID, ' }\n');
fprintf(fileID, '}\n');

% probes
file_probes_output = 'output.out';
fprintf(fileID, 'nearfield %s %s\n', file_probes_input, file_probes_output);

disp('objective configuration created');
%% simulated fields
call_exe('std_config');

%% load E2(x1)
data = load(file_probes_output);
make_complex = @(x, y) x + 1j * y;

E2 = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H2 = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E2 = E2 / ref_signal_fft;
H2 = H2 / ref_signal_fft;

%% E1(x2) * G2
% write to probes file
file_probes_input = 'probes.in';
fileID = fopen(file_probes_input, 'w');
fprintf(fileID, '1\n');
fprintf(fileID, '%e %e %e %e\n', x2, fp);
fclose(fileID);

% write to dipoles file
ctype = [1, 2, 4];
file_dipoles_input = 'dipoles.in';
fileID = fopen(file_dipoles_input, 'w');
fprintf(fileID, '3\n');
for i = 1 : 3
    fprintf(fileID, '%e %e %e %e %e %e %d\n', x1, abs(G1(i)), fp, -angle(G1(i)) / (2 * pi * fp), ctype(i));
end
fclose(fileID);

%% simulated fields
call_exe('std_config');

%% load E1(x2)
data = load(file_probes_output);
make_complex = @(x, y) x + 1j * y;

E1 = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H1 = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E1 = E1 / ref_signal_fft;
H1 = H1 / ref_signal_fft;
%% Calculation

inner_prod1 = sum(G1 .* E2);
inner_prod2 = sum(G2 .* E1);
disp(inner_prod1)
disp(inner_prod2)

