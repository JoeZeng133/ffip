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
ft = c0 / (500e-9); %frequency to extract fields at
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
inhom.len = [80e-9, 80e-9, 12e-9];
inhom.center = [center(1), center(2), 45e-9];
inhom.radius = 30e-9;
inhom.shape = make_disk([0 0], inhom.radius, 'num_elements', 40);
[inhom.rho, inhom.x, inhom.y, inhom.z, inhom.p1, inhom.p2, inhom.dim, inhom.dx] = ...
    make_inhom(inhom.shape, inhom.len(3), -inhom.len(1:2)/2, inhom.len(1:2)/2, 'dx', 4e-9, 'center', inhom.center, 'float', true);
inhom.position = inhom.p1;
inhom.type = 'inhom';
inhom.filename = 'inhom.in';
inhom.medium_idx1 = 1;
inhom.medium_idx2 = 0;
inhom.num = numel(inhom.rho);

[~, er1] = Au(2 * pi * ft);
[~, er2] = Air(2 * pi * ft);

% target probes
target.p1 = [0 0 60e-9];
target.p2 = [p2(1) p2(2) 60e-9];
[target.x, target.y, target.z, target.dim, target.dx] = make_interp_region(target.p1, target.p2, 'dim', [11 11 1]);
target.freq = ft * ones(target.dim);
target.num = numel(target.x);

% nearfield probes at target positions
nf.x = target.x(:);
nf.y = target.y(:);
nf.z = target.z(:);
nf.freq = ft * ones(target.dim);
nf.input_file = 'nf.in';
nf.output_file = 'nf.out';

%% write configuration
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
call_exe('std_config');

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

disp('objectiev fields calcualted');
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
%% 
rho = ones(inhom.dim(1:2));                          %start with rho=1
tol = 1e-6;
step = 0.1; 
a = 0.9;
max_itr = 2;
itr = 0;
max_size = [1 1000];
inhom_forward = inhom;

% logging data
obj_func = zeros(max_size);
exp_obj_func = zeros(max_size);
sens = cell(max_size);
Se_list = cell(max_size);
rho_list = cell(max_size);

%% optimization
while (itr < max_itr)
    itr = itr + 1;
    fprintf('################# Iteration %d\n', itr);
    
    if (mod(itr - 1, 1) == 0)
        figure
        surf(inhom.x(:, :, 1)/1e-9, inhom.y(:, :, 1)/1e-9, rho)
        colorbar
        shading flat
        caxis([0 1])
        zlim([0 1])
        title(['Iteration ', num2str(itr)])
        xlabel('x')
        ylabel('y')
        getframe;
    end
    
    % forward geometry created
    rho_list{itr} = rho;
    inhom_forward.rho = rho .* ones(inhom.dim);

    % nf probes (inhom + target)
    forward_output_file = 'nf_forward.out';
    nf_forward.x = [inhom.x(:) ; target.x(:)];
    nf_forward.y = [inhom.y(:) ; target.y(:)];
    nf_forward.z = [inhom.z(:) ; target.z(:)];
    nf_forward.freq = ft * ones([inhom.num + target.num, 1]);
    nf_forward.input_file = 'nf_forward.in';
    nf_forward.output_file = forward_output_file;
    
    geometry = {inhom_forward, fept};
    source = {plane_wave};
    gen_config(basic, medium, geometry, source, 'nearfield', nf_forward, 'step_output', 0);
    disp('forward simulation configuration complete');
    % forward simulation results
    call_exe('std_config')
    data = load(forward_output_file);

    E_forward = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
    H_forward = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
    E_forward = E_forward / ref_plane_fft;
    H_forward = H_forward / ref_plane_fft;

    E_forward_inhom = E_forward(1 : inhom.num, :);
    H_forward_inhom = H_forward(1 : inhom.num, :);
    E_forward_target = E_forward(inhom.num + 1 : end, :);
    H_forward_target = H_forward(inhom.num + 1 : end, :);
    
    disp('forward calculation complete')

%     target_func = 0.5 * sum(abs(E_forward_probes(:, 3) - E_target_probes(:, 3)).^2, 1); %only Ez
    obj_func(itr) = 0.5 * sum(sum(abs(E_forward_target - E_target).^2, 2), 1);
%     obj_func(itr) = 0.5 * sum(sum(abs(H_forward_target - H_target_probes).^2, 2), 1);
    fprintf('Current objective function is evaluated at : %e\n', obj_func(itr));
    
    % adjoint simulation configuration
    ref_dipole = ricker(t, fs, delay);
%     ref_dipole = sin(2 * pi * fs * t);
    ref_dipole_fft = sum(ref_dipole.* exp(-1j * 2 * pi * ft * t), 2);
    
    Ex = 1;
    Ey = 2;
    Ez = 4;
    Hx = 7 - Ex;
    Hy = 7 - Ey;
    Hz = 7 - Ez;
    dipole.data = conj(E_forward_target - E_target);
    dipole.template = ones(size(dipole.data));
    dipole.x = repmat(target.x(:), 1, 3);
    dipole.y = repmat(target.y(:), 1, 3);
    dipole.z = repmat(target.z(:), 1, 3);
    dipole.func_type = 'r' * dipole.template;
    dipole.delay = delay - angle(dipole.data) / (2 * pi * ft);
    dipole.amp = abs(dipole.data);
    dipole.fp = fs * dipole.template;
    dipole.ctype = [Ex Ey Ez] .* dipole.template;
    dipole.type = 'dipole';
    dipole.filename = 'adjoint_dipoles.in';
    
    geometry = {inhom_forward, fept};
    source = {dipole};
    adjoint_output_file = 'nf_adjoint.out';
    nf_forward.output_file = adjoint_output_file;
    gen_config(basic, medium, geometry, source, 'nearfield', nf_forward, 'step_output', 0);
    
    % adjoint fields
    call_exe('std_config');
    data = load(adjoint_output_file);
    E_adjoint = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
    H_adjoint = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
    E_adjoint = E_adjoint / ref_dipole_fft;
    H_adjoint = H_adjoint / ref_dipole_fft;

    E_adjoint_inhom = E_adjoint(1 : inhom.num, :);
    H_adjoint_inhom = H_adjoint(1 : inhom.num, :);
    E_adjoint_target = E_adjoint(inhom.num + 1 : end, :);
    H_adjoint_target = H_adjoint(inhom.num + 1 : end, :);

    disp('adjoint simulation completed');
    % Sensitivity Calculation    
    [w_X, w_Y, w_Z] = ndgrid(...
        make_2segment(0.5, 1, inhom.dim(1)), ...
        make_2segment(0.5, 1, inhom.dim(2)), ...
        make_2segment(0.5, 1, inhom.dim(3)));
    w = w_X .* w_Y .* w_Z;
    w = w(:) * prod(inhom.dx);
    
    Se = 1j * (2 * pi * ft) * sum(E_forward_inhom .* E_adjoint_inhom, 2);
    de_drho = e0 * (er1 - er2);
    A = real(Se .* de_drho * prod(inhom.dx));
    A = reshape(A, inhom.dim);
    A = sum(A, 3);

    Se_list{itr} = Se;
    sens{itr} = A;
    
    disp('sensitivity calculation complete')
    % Optimization
    lb = max(-rho(:), -step);
    ub = min(1 - rho(:), step);
    delta_rho = linprog(A(:), [], [], [], [], lb, ub);
    exp_obj_func(itr) = obj_func(itr) + sum(delta_rho .* A(:));
    fprintf('The new objective function is expected to be : %.2e\n', exp_obj_func(itr));

    % updates
    rho = reshape(delta_rho, size(rho)) + rho;
    step = a * step;
    disp('updated');
end




