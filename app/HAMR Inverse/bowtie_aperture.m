close all
clear
clc
fclose all;

%% build geometry
a = 30e-9;
b = 40e-9;
c = 50e-9;
d = 50e-9;

tri1 = polyshape([0, c, -c], [b, b + d, b + d]);
tri2 = polyshape([0, c, -c], [-b, -b - d, -b - d]);
rect = make_rect([0 0], [a, (b + d) * 2]);

tmp = union(tri1, tri2);
au_shape = union(tmp, rect);

bbox_p1 = [-150e-9, -150e-9];
bbox_p2 = [150e-9, 150e-9];
bbox = make_rect((bbox_p1 + bbox_p2) / 2, (bbox_p2 - bbox_p1));

figure
plot(scale(au_shape, 1e9)), hold on
plot(scale(bbox, 1e9)), hold off
xlabel('x [nm]')
ylabel('y [nm]')
axis equal

%% build configuration
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

er_bg = 1;
ur_bg = 1;
PML_d = 6;
Sc = 0.5;
dx = 2e-9;
dt = Sc * dx / c0;
dim = [150, 150, 50];
step = 10000;
tf_layer = 2;
sf_layer = 1;
projector_padding = ceil((tf_layer + 1) / 2);
p1 = [0 0 0];
p2 = dim * dx;
center = dim * dx / 2;

fs = c0 / 800e-9;
delay = 1 / fs;
ft = c0 / 800e-9;

% FePt material
FePt = sing_freq((2.9 - 1.5j)^2, ft);

% gold aperture
[au_rho, au_x, au_y, au_z, au_p1, au_p2, au_cell] = make_inhom(au_shape, 60e-9, bbox_p1, bbox_p2, dx);
au_center = [center(1); center(2); 35e-9];
au_geom = struct('type', 'inhom', 'medium_idx1', 0, 'medium_idx2', 1, 'filename', 'au_geom.in',...
    'dim', size(au_x), 'position', au_p1 + au_center, 'dx', au_cell, 'rho', au_rho * 1);

%FePt layer
fe_center = [center(1) center(2) 75e-9];
fe_len = [dim(1) * dx, dim(2) * dx, 10e-9];
fe_p1 = fe_center - fe_len / 2;
fe_p2 = fe_center + fe_len / 2;
fe_geom = struct('type', 'box', 'medium_idx', 2, 'lower_position', fe_p1, 'upper_position', fe_p2);

% gold fields
[au_nf.x, au_nf.y, au_nf.z] = make_interp_region([0, 0, au_center(3)], [p2(1), p2(2), au_center(3)], dx * 2);
au_nf.freq = ft * ones(size(au_nf.x));
num_au_nf = numel(au_nf.x);

% fept fields
[fe.x, fe.y, fe.z] = make_interp_region([0, 0, fe_center(3)], [p2(1), p2(2), fe_center(3)], dx * 2);
fe.freq = ft * ones(size(fe.x));
num_fe_nf = numel(fe.x);

%% write configuration
basic.er_bg = er_bg;
basic.ur_bg = ur_bg;
basic.PML_d = PML_d;
basic.dx = dx;
basic.dt = dt;
basic.dim = dim;
basic.step = step;
basic.tf_layer = tf_layer;
basic.sf_layer = sf_layer;

medium{1} = Air();
medium{2} = Au();
medium{3} = FePt;

geometry = {au_geom, fe_geom};

source{1} = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos',...
    dim(3) + projector_padding, 'func_type', 'r', 'fp', fs, 'delay', delay, 'ref_pos', 0);

nf = struct('input_file', 'nf.in', 'output_file', 'output.out', 'x', [au_nf.x(:) ; fe.x(:)], ...
    'y', [au_nf.y(:) ; fe.y(:)], 'z', [au_nf.z(:) ; fe.z(:)], 'freq', [au_nf.freq(:) ; fe.freq(:)]);

gen_config(basic, medium, geometry, source, 'nearfield', nf, 'step_output', 1);


%%
ref_plane_signal = load('reference.out');
ref_plane_signal = reshape(ref_plane_signal, 1, []);
t = (0:numel(ref_plane_signal) - 1) * dt;
ref_plane_signal_fft = sum(ref_plane_signal .* exp(-1j * 2 * pi * ft(:) * t), 2);

data = load(nf.output_file);
make_complex = @(x, y) x + 1j * y;

E = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E = E / ref_plane_signal_fft;
H = H / ref_plane_signal_fft;

E_au = reshape(E(1:num_au_nf, :), [size(au_nf.x), 3]);
H_au = reshape(H(1:num_au_nf, :), [size(au_nf.x), 3]);
E_au_abs = sqrt(abs(E_au(:, :, 1)).^2 + abs(E_au(:, :, 2)).^2 + abs(E_au(:, :, 3)).^2);

figure
surf(au_nf.x / 1e-9, au_nf.y / 1e-9, abs(E_au(:, :, 1)))
xlabel('x')
ylabel('y')
shading flat
colormap jet
colorbar
title('|E| in Au layer')

E_fe = reshape(E(num_au_nf+1:end, :), [size(fe.x), 3]);
H_fe = reshape(H(num_au_nf+1:end, :), [size(fe.x), 3]);
E_fe_abs = sqrt(abs(E_fe(:, :, 1)).^2 + abs(E_fe(:, :, 2)).^2 + abs(E_fe(:, :, 3)).^2);

figure
surf(fe.x / 1e-9, fe.y / 1e-9, E_fe_abs)
xlabel('x')
ylabel('y')
shading flat
colormap jet
colorbar
title('|E| in FePt layer')
