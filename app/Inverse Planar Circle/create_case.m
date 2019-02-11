clear
clc
close all
fclose all;

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
fs = c0 / (500e-9); %center frequency of the pulse
ft = c0 / (800e-9); %frequency to extract fields at
delay = 1 / fs;
fprintf('Equal to %.2f times of full pulse width\n', time_step * dt * fs);

center = dim / 2 * dx;
p2 = dim * dx;

% FePt 
fept.len = [100e-9, 100e-9, 10e-9];
fept.center = [center(1), center(2), 60e-9];
fept.lower_position = fept.center - fept.len / 2;
fept.upper_position = fept.center + fept.len / 2;
fept.type = 'box';
fept.medium_idx = 2;

% disk antenna
inhom.len = [70e-9, 70e-9, 10e-9];
inhom.center = [center(1), center(2), 45e-9];
inhom.radius = 30e-9;
inhom.shape = make_disk([0 0], inhom.radius, 'num_elements', 40);
[inhom.rho, inhom.x, inhom.y, inhom.z, inhom.p1, inhom.p2, inhom.dim, inhom.dx] = ...
    make_inhom(inhom.shape, inhom.len(3), -inhom.len(1:2)/2, inhom.len(1:2)/2, 'dx', 4e-9, 'center', inhom.center, 'float', true);
inhom.position = inhom.p1;
inhom.type = 'inhom';
inhom.filename = 'disk.in';
inhom.medium_idx1 = 1;
inhom.medium_idx2 = 0;
inhom.num = numel(inhom.rho);

% target probes
target.p1 = [0 0 60e-9];
target.p2 = [p2(1) p2(2) 60e-9];
[target.x, target.y, target.z, target.dim, target.dx] = make_interp_region(target.p1, target.p2, 'dim', [11 11 1]);
target.freq = ft * ones(target.dim);
target.num = numel(target.x);

% nearfield probes at target positions
nf.x = [target.x(:)];
nf.y = [target.y(:)];
nf.z = [target.z(:)];
nf.freq = ft * ones(target.dim);
nf.input_file = 'nf.in';
nf.output_file = 'nf.out';

%%
basic.er_bg = er_bg;
basic.ur_bg = ur_bg;
basic.PML_d = PML_d;
basic.dx = dx;
basic.dt = dt;
basic.dim = dim;
basic.step = time_step;
basic.tf_layer = tf_layer;
basic.sf_layer = sf_layer;

medium = {Air(), Au(), sing_freq((2.9 - 1.5j)^2, ft)}; %Air, Au, FePt
geometry = {inhom, fept};
plane_wave = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos',...
    dim(3) + projector_padding, 'func_type', 'r', 'fp', fs, 'delay', delay, 'ref_pos', 0); 
source = {plane_wave};
gen_config(basic, medium, geometry, source, 'nearfield', nf, 'step_output', 1);
%% simulated fields
% call_exe('std_config');

%% load data
ref_plane = load('reference.out');
ref_plane = reshape(ref_plane, 1, []);
t = (0:numel(ref_plane) - 1) * dt;
ref_plane_fft = sum(ref_plane .* exp(-1j * 2 * pi * ft * t), 2);

data = load(nf.output_file);
make_complex = @(x, y) x + 1j * y;

E_target = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H_target = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];

E_target = E_target / ref_plane_fft;
H_target = H_target / ref_plane_fft;

save('case_configuration');
disp('objective fields saved');

%% visualization
figure
plot(inhom.shape), hold on
plot(make_rect([0 0], inhom.len(1:2))), hold off
xlabel('x')
ylabel('y')

figure
H = sqrt(sum(abs(H_target).^2, 2));
H = reshape(H, target.dim);
surf(target.x, target.y, H)
shading flat
xlabel('x')
ylabel('y')
title('|H|')

figure
E = sqrt(sum(abs(E_target).^2, 2));
E = reshape(E, target.dim);
surf(target.x, target.y, E)
shading flat
xlabel('x')
ylabel('y')
title('|E|')
