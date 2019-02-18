%% plane wave test
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

er_bg = 1;
ur_bg = 1;
PML_d = 6;
Sc = 0.5;
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [30, 30, 30];
step = 600;
tf_layer = 0;
sf_layer = 0;
projector_padding = ceil((tf_layer + 1) / 2);
center = dim * dx / 2;
p2 = dim * dx;

fs = c0 / (30 * dx);
delay = 1 / fs;

lambda = 30 * dx;
ft = c0 / (30 * dx);
a = 10 * dx;
%% generate configuration
basic.er_bg = er_bg;
basic.ur_bg = ur_bg;
basic.PML_d = PML_d;
basic.dx = dx;
basic.dt = dt;
basic.dim = dim;
basic.step = step;
basic.tf_layer = tf_layer;
basic.sf_layer = sf_layer;

medium = {fic1()};

sphere = struct('type', 'sphere', 'medium_idx', 0, 'radius', a, 'position', [12 12 15] * dx);
box = struct('type', 'box', 'lower_position', [-1 -1 5 * dx], 'upper_position', [1 1 10 * dx], 'medium_idx', 0);
geometry = {sphere};

source = {struct('type', 'plane2', 'func_type', 'r', 'fp', fs, 'delay', delay, 'pos', dx)};

nf.p1 = [center(1) 0 0];
nf.p2 = [center(1) p2(2) p2(3)];
[nf.x, nf.y, nf.z, nf.dim, nf.dx] = make_interp_region(nf.p1, nf.p2, 'dx', dx);
nf.freq = ft * ones(nf.dim);
nf.input_file = 'nf.in';
nf.output_file = 'output.out';



gen_config(basic, medium, geometry, source, 'step_output', 1, 'num_proc', 2);

disp('config.in created');

%%
file = fopen('snapshot.out', 'r');
dim = fscanf(file, '%e %e %e', [1 3]);
data = fscanf(file, '%e');
fclose(file);


data = reshape(data, dim(2), dim(3), []);
norm = max(abs(data(:)));
interval = 0.1;

tic;

s = toc;
for i = 1 : size(data, 3)
    surf(data(:, :, i) / norm);
    xlabel('z')
    ylabel('y')
    zlim([-1 1])
    
    pause(interval)
    s = toc;
    getframe;
end

%%
surf(data(:, :, 100) / norm);
xlabel('z')
ylabel('y')

