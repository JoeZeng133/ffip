%% mie theory test, assuming vaccum background 
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
dim = [40, 40, 40];
step = 1000;
tf_layer = 2;
sf_layer = 1;
projector_padding = ceil((tf_layer + 1) / 2);

fs = c0 / (30 * dx);
delay = 1 / fs;

lambda = 30 * dx;
ft = c0 / (30 * dx);
a = 10 * dx;

center = dim * dx / 2;
rho = a + linspace(2 * dx, 10 * dx, 10);
phi = linspace(0, 2 * pi, 10);
th = linspace(0, pi, 10);
[xc, yc, zc] = sph2cart(phi, pi / 2 - th, rho);
[Xc, Yc, Zc, Ft] = ndgrid(xc, yc, zc, ft);

Omega = 2 * pi * Ft;
K = Omega / c0;
er_func = @fic1;

[~, er] = er_func(2 * pi * ft);
ns = sqrt(conj(er));
nm = 1;

% inhomogeneous region
inhom.center = center;
inhom.length = [30 * dx, 30 * dx, 30 * dx];
inhom.position = inhom.center - inhom.length / 2;
inhom.dim = round(inhom.length / dx) + 1;
inhom.dx = inhom.length ./ (inhom.dim - 1);
[inhom.x, inhom.y, inhom.z] = ndgrid(...
    inhom.position(1) + (0:inhom.dim(1) - 1) * inhom.dx(1), ...
    inhom.position(2) + (0:inhom.dim(2) - 1) * inhom.dx(2), ...
    inhom.position(3) + (0:inhom.dim(3) - 1) * inhom.dx(3));

inhom.rho = (inhom.x - center(1)).^2 + (inhom.y - center(2)).^2 + (inhom.z - center(3)).^2 <= a^2;
inhom.rho = inhom.rho * 1;
inhom.type = 'inhom';
inhom.medium_idx1 = 0;
inhom.medium_idx2 = 1;
inhom.filename = 'inhom.in';

% blobs
[blob.x, blob.y, blob.z] = ndgrid((0:dim(1)) * dx, (0:dim(2)) * dx, (0:dim(3)) * dx);
blob.amp = ((blob.x - center(1)).^2 + (blob.y - center(2)).^2 + (blob.z - center(3)).^2) <= a^2;
blob.x = blob.x(blob.amp);
blob.y = blob.y(blob.amp);
blob.z = blob.z(blob.amp);
blob.amp = blob.amp(blob.amp) * 1;
blob.medium_idx = zeros(size(blob.amp));
blob.input_file = 'blob.in';
blob.type = 'blob';


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

medium{1} = fic1();
medium{2} = Air();

% geometry{1} = struct('type', 'sphere', 'medium_idx', 0, 'radius', a, 'position', dim * dx / 2);
geometry{1} = inhom;
% geometry{1} = blob;

source{1} = struct('type', 'plane', 'dim_neg', projector_padding, 'dim_pos', ...
    dim(3) + projector_padding, 'func_type', 'r', 'fp', fs, 'delay', delay, 'ref_pos', dim(3) / 2 * dx);


nf.x = Xc + center(1);
nf.y = Yc + center(2);
nf.z = Zc + center(3);
nf.freq = Ft;
nf.input_file = 'nf.in';
nf.output_file = 'output.out';

gen_config(basic, medium, geometry, source, 'nearfield', nf, 'step_output', 1, 'num_proc', 2);

disp('config.in created');


%% theoretical fields
[E, H] = calcmie_nf( a, ns, nm, lambda, Xc(:), Yc(:), Zc(:), 'TotalField', true );

%% numerical fields
call_exe('std_config')
data = load('output.out');
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
    figure
    subplot(2, 2, 1)
    plot(real(E(:, i)), real(En(:, i)), '.'), hold on
    plot(real(E(:, i)), real(E(:, i))), hold off
    xlabel('Analytical')
    axis equal
    title(num2str(i))
    
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


