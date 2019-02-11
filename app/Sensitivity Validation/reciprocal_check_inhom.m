%% blob check
clear
clc
close all
fclose all;

%% basic configuration
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
dim = [30, 30, 30];
step = 10000;
tf_layer = 2;
sf_layer = 1;
projector_padding = ceil((tf_layer + 1) / 2);
center = dim * dx / 2;
p1 = [0 0 0];
p2 = dim * dx;

% excitation frequency
fs = c0 / 800e-9;
delay = 0;

% frequency of interests
lambda = 800e-9;
ft = c0 / lambda;

[~, er] = Au(2 * pi * ft);

inhom.p1 = [20e-9, 20e-9, 20e-9];
inhom.p2 = [40e-9, 40e-9, 40e-9];
[inhom.x, inhom.y, inhom.z, inhom.dim, inhom.dx] = make_interp_region(inhom.p1, inhom.p2, 2 * dx);
inhom.rho = 1 * ones(inhom.dim);
rho_copy = inhom.rho;
inhom.position = inhom.p1;
inhom.medium_idx1 = 1;
inhom.medium_idx2 = 0;
inhom.filename = 'inhom.in';
inhom.type = 'inhom';
num_inhom = numel(inhom.x);

nf.p1 = [10e-9, 10e-9, 50e-9];
nf.p2 = [50e-9, 50e-9, 50e-9];
[nf.x, nf.y, nf.z, nf.dim] = make_interp_region(nf.p1, nf.p2, dx);
num_ob = numel(nf.x);
nf.x = [nf.x(:) ; inhom.x(:)];
nf.y = [nf.y(:) ; inhom.y(:)];
nf.z = [nf.z(:) ; inhom.z(:)];
nf.freq = ft * ones(size(nf.x));
nf.input_file = 'nf.in';
nf.output_file = 'output.out';
%% original configuration
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

sphere = struct('type' ,'sphere', 'radius', 15e-9, 'position', center, 'medium_idx', 1);
blob_sphere = struct('type', 'sphere', 'radius', 2e-9, 'position', center, 'medium_idx', 0);
box = struct('type', 'box', 'lower_position', inhom.p1, 'upper_position', inhom.p2, 'medium_idx', 1);
layer = struct('type', 'box', 'lower_position', [p1(1) p1(2) inhom.p1(3)], 'upper_position', [p2(1) p2(2) inhom.p2(3)], 'medium_idx', 1);
blob_region = struct('type', 'box', 'lower_position', inhom.p1, 'upper_position', inhom.p2, 'medium_idx', 1);

geometry = {inhom};
plane_wave = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos', ...
    dim(3) + projector_padding, 'func_type', 's', 'fp', fs, 'delay', delay, 'ref_pos', dim(3) / 2 * dx);

source = {plane_wave};

gen_config(basic, medium, geometry, source, 'nearfield', nf);
disp('config.in created');

%% original fields
call_exe('std_config');
data = load(nf.output_file);
make_complex = @(x, y) x + 1j * y;
ref_plane = load('reference.out');
ref_plane = reshape(ref_plane, 1, []);
t = (0:numel(ref_plane) - 1) * dt;
ref_plane_fft = sum(ref_plane .* exp(-1j * 2 * pi * nf.freq(:) * t), 2);

Eo = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
Ho = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
Eo = Eo ./ ref_plane_fft;
Ho = Ho ./ ref_plane_fft;

disp('Original fields obtained'); 
%% perturbed configuration
[px, py, pz] = ndgrid(make_2segment(0, 1, inhom.dim(1)), make_2segment(0, 1, inhom.dim(2)), make_2segment(0, 1, inhom.dim(3)));
weight = px .* py .* pz;
pert = -0.1 * weight .* rand(inhom.dim);
inhom.rho = rho_copy + pert;
geometry = {inhom};
source = {plane_wave};
gen_config(basic, medium, geometry, source, 'nearfield', nf);
disp('config.in created');

%% perturbed fields
call_exe('std_config');
data = load(nf.output_file);
Ep = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
Hp = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
Ep = Ep ./ ref_plane_fft;
Hp = Hp ./ ref_plane_fft;
disp('perturbed fields obtained')

%% adjoint configuration
inhom.rho = rho_copy;

Ex = 1;
Ey = 2;
Ez = 4;

dipole.data = 1j * (2 * pi * ft) * e0 * (er - 1) * pert(:) .* weight(:) .* Eo(num_ob + 1 : end, :) * prod(inhom.dx);
dipole.template = ones(size(dipole.data));
dipole.x = repmat(inhom.x(:), 1, 3);
dipole.y = repmat(inhom.y(:), 1, 3);
dipole.z = repmat(inhom.z(:), 1, 3);
dipole.func_type = 's' * dipole.template;
dipole.delay = - angle(dipole.data) / (2 * pi * ft);
dipole.amp = abs(dipole.data);
dipole.fp = fs * dipole.template;
dipole.ctype = [Ex Ey Ez] .* dipole.template;
dipole.type = 'dipole';
dipole.filename = 'dipole.in';

ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
ref_dipole = sin(2 * pi * fs * t);
ref_dipole_fft = sum(ref_dipole.* exp(-1j * 2 * pi * ft * t), 2);

geometry = {inhom};
source = {dipole};
gen_config(basic, medium, geometry, source, 'nearfield', nf);
disp('config.in created');

%% adjoint field
call_exe('std_config');
data = load(nf.output_file);
Ea = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
Ha = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
Ea = Ea ./ ref_dipole_fft;
Ha = Ha ./ ref_dipole_fft;
disp('adjoint fields obtained')

%%
Ediff = Ep - Eo;

Ediff_ob = Ediff(1:num_ob, :);
Ea_ob = Ea(1:num_ob, :);

Ediff_blob = Ediff(num_ob + 1:end, :);
Ea_blob = Ea(num_ob + 1:end, :);

figure
plot(imag(Ediff_ob(:)), imag(Ea_ob(:)), '.'), hold on
plot(imag(Ediff_ob(:)), imag(Ediff_ob(:))), hold off
title('imag E')

figure
plot(real(Ediff_ob(:)), real(Ea_ob(:)), '.'), hold on
plot(real(Ediff_ob(:)), real(Ediff_ob(:))), hold off
title('real E')
