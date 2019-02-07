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
fs = c0 / 500e-9;
delay = 0;

% frequency of interests
lambda = 500e-9;
ft = c0 / lambda;

[~, er] = fic1(2 * pi * ft);

blob.p1 = [20e-9, 20e-9, 20e-9];
blob.p2 = [40e-9, 40e-9, 40e-9];
[blob.x, blob.y, blob.z, blob.dim] = make_interp_region(blob.p1, blob.p2, dx);
blob.amp = -0.01 * ones(blob.dim);
blob.medium_idx = 1 * ones(blob.dim);
blob.input_file = 'blob.in';
blob.type = 'blob';
num_blob = numel(blob.x);

nf.p1 = [10e-9, 10e-9, 30e-9];
nf.p2 = [50e-9, 50e-9, 40e-9];
[nf.x, nf.y, nf.z, nf.dim] = make_interp_region(nf.p1, nf.p2, dx);
num_ob = numel(nf.x);
nf.x = [nf.x(:) ; blob.x(:)];
nf.y = [nf.y(:) ; blob.y(:)];
nf.z = [nf.z(:) ; blob.z(:)];
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
medium{2} = fic1();

gold_box = struct('type', 'box', 'lower_position', blob.p1, 'upper_position', blob.p2, 'medium_idx', 1);
gold_layer = struct('type', 'box', 'lower_position', [p1(1) p1(2) blob.p1(3)], 'upper_position', [p2(1) p2(2) blob.p2(3)], 'medium_idx', 1);
blob_region = struct('type', 'box', 'lower_position', blob.p1, 'upper_position', blob.p2, 'medium_idx', 1);

geometry = {blob_region};
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

disp('Fields Extracted'); 
%% perturbed configuration
geometry = {blob_region, blob};
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

%% adjoint configuration
Ex = 1;
Ey = 2;
Ez = 4;

dipole.data = 1j * (2 * pi * ft) * e0 * er * blob.amp(:) .* Eo(num_ob + 1 : end, :) * dx^3;
dipole.template = ones(size(dipole.data));
dipole.x = repmat(blob.x(:), 1, 3);
dipole.y = repmat(blob.y(:), 1, 3);
dipole.z = repmat(blob.z(:), 1, 3);
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

geometry = {blob_region};
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
