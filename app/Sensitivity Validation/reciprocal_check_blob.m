%% mie theory test, assuming vaccum background 
clear
clc
close all
fclose all;

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

fs = c0 / 500e-9;
delay = 1 / fs;

lambda = 500e-9;
ft = c0 / lambda;
a = 40e-9;

center = dim * dx / 2;
Omega = 2 * pi * Ft;
K = Omega / c0;
[~, er] = Au(2 * pi * ft);

ns = sqrt(conj(er));
nm = 1;
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

medium{1} = Air();
medium{2} = Au();

geometry = {};

source{1} = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos', ...
    dim(3) + projector_padding, 'func_type', 'r', 'fp', fs, 'delay', delay, 'ref_pos', dim(3) / 2 * dx);

nf.x = Xc + center(1);
nf.y = Yc + center(2);
nf.z = Zc + center(3);
nf.freq = Ft;
nf.input_file = 'nf.in';
nf.output_file = 'output.out';

gen_config(basic, medium, geometry, source, 'nearfield', nf);

disp('config.in created');

%% numerical fields
data = load(nf.output_file);
make_complex = @(x, y) x + 1j * y;
ref_signal = load('reference.out');
ref_signal = reshape(ref_signal, 1, []);
t = (0:numel(ref_signal) - 1) * dt;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

En = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
Hn = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
En = En ./ ref;
Hn = Hn ./ ref;

disp('Numerical Fields Extracted');


%% corellation plots
for i = 1 : 3
    figure(i + 1)
    subplot(2, 2, 1)
    plot(real(E(:, i)), real(En(:, i)), '.'), hold on
    plot(real(E(:, i)), real(E(:, i))), hold off
    xlabel('Analytical')
    axis equal
    
    subplot(2, 2, 2)
    plot(imag(E(:, i)), -imag(En(:, i)), '.'), hold on
    plot(imag(E(:, i)), imag(E(:, i))), hold off
    xlabel('Analytical')
    axis equal

    subplot(2, 2, 3)
    plot(real(H(:, i)), real(Hn(:, i)), '.'), hold on
    plot(real(H(:, i)), real(H(:, i))), hold off
    xlabel('Analytical')
    axis equal
    
    subplot(2, 2, 4)
    plot(imag(H(:, i)), -imag(Hn(:, i)), '.'), hold on
    plot(imag(H(:, i)), imag(H(:, i))), hold off
    xlabel('Analytical')
    axis equal
end

