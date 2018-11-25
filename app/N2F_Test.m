%% far-field dipole test
clear
clc
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);

Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [10, 10, 10];
step = 600;

Np = 20;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
ricker = @(t, d, fp) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
t = (0:step) * dt;
ref_signal = ricker(t, 1/fp, fp);

rho = (1e4 + 0.33) * (Np * dx);
phi = pi;
th = linspace(0.25 * pi, 0.75 * pi, 100);
ft = fp;

[Ft, Th, Phi, Rho] = ndgrid(ft, th, phi, rho);

pos = [Ft(:), Th(:), Phi(:), Rho(:)];

fileID = fopen('request.in', 'w');
fprintf(fileID, '%d\n', size(pos, 1));
fprintf(fileID, '%e %e %e %e\n', pos');
fclose(fileID);


%% theoretical fields
Omega = 2 * pi * Ft;

er =  1;
e = e0 * er;
ur = 1;
u = u0;

c = 1 ./ sqrt(e * u);
eta = sqrt(u ./ e);
K = Omega ./ c;

Hphi_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * K ./ Rho + 1 ./ Rho.^2) .* sin(Th);
Er_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (2 * eta ./ Rho.^2 + 2 ./ (1j * Omega .* e .* Rho.^3)) .* cos(Th);
Eth_p = 1 / (4 * pi) * exp(-1j * K .* Rho) .* (1j * Omega .* u ./ Rho + 1 ./ (1j * Omega .* e .* Rho.^3) + eta ./ Rho.^2) .* sin(Th);


%% simulated fields
data = load('data.out');
make_complex = @(x, y) x + 1j * y;
ref = sum(ref_signal .* exp(-1j * 2 * pi * Ft(:) * t), 2);

Eth = make_complex(data(:, 1), data(:, 2)) ./ ref;
Ephi = make_complex(data(:, 3), data(:, 4)) ./ ref;
Hth = make_complex(data(:, 5), data(:, 6)) ./ ref;
Hphi = make_complex(data(:, 7), data(:, 8)) ./ ref;

%% comparisons
fprintf('Ephi over Eth = %e \n', max(abs(Ephi)) / max(abs(Eth)))
fprintf('Hth over Hphi = %e \n', max(abs(Hth)) / max(abs(Hphi)))
figure(2)
subplot(1, 2, 1)
plot(real(Hphi(:)), 'r'), hold on
plot(real(Hphi_p(:)), 'b'), hold off
title('Re H_\phi')

subplot(1, 2, 2)
plot(imag(Hphi(:)), 'r'), hold on
plot(imag(Hphi_p(:)), 'b'), hold off
title('Im H_\phi')

figure(3)
subplot(1, 2, 1)
plot(real(Eth(:)), 'r'), hold on
plot(real(Eth_p(:)), 'b'), hold off
title('Re E_\theta')

subplot(1, 2, 2)
plot(imag(Eth(:)), 'r'), hold on
plot(imag(Eth_p(:)), 'b'), hold off
title('Im E_\theta')

%% correlation comparisons
% figure(4)
% subplot(1, 2, 1)
% plot(real(Eth(:)), real(Eth_p(:)), '*')
% axis equal
% title('real E_\theta')
% 
% subplot(1, 2, 2)
% plot(imag(Eth(:)), imag(Eth_p(:)), '*')
% axis equal
% title('imag E_\theta')
% 
% figure(5)
% subplot(1, 2, 1)
% plot(real(Hphi(:)), real(Hphi_p(:)), '*')
% axis equal
% title('real H_\phi')
% 
% subplot(1, 2, 2)
% plot(imag(Hphi(:)), imag(Hphi_p(:)), '*')
% axis equal
% title('imag H_\phi')

%% manually read fields at each time and convert to N, L to far fields
clc

data_raw = cell([6 1]);
norm_vec = zeros([6 3]);
int_weights = cell([6 1]);
pos = cell([6 1]);
E = cell([6 1]);
H = cell([6 1]);
ref_signal2 = ricker((1:step-1) * dt, 1/fp, fp);
ref_f = fft(ref_signal2);

for i = 1 : 6
    filename = ['face', num2str(i - 1), '.out'];
    fs = fopen(filename);
    p1 = fscanf(fs, '%e %e %e', [1 3]);
    p2 = fscanf(fs, '%e %e %e', [1 3]);
    norm_vec(i, :) = fscanf(fs, '%e %e %e', [1 3]);
    face_dim = (p2 - p1) / 2 + 1;
    
    % get integration weight for each point
    int = cell([3 1]);
    for j =  1 : 3
        int{j} = ones([face_dim(j), 1]);
        if(face_dim(j) ~= 1)
            int{j}(1) = 0.5;
            int{j}(end) = 0.5;
        end    
    end
    [intx, inty, intz] = ndgrid(int{1}, int{2}, int{3});
    int_weights{i} = intx .* inty .* intz;
    int_weights{i} = int_weights{i}(:);
    
    % get position - center for each point
    [X, Y, Z] = ndgrid(p1(1):2:p2(1), p1(2):2:p2(2), p1(3):2:p2(3));
    pos{i} = ([X(:), Y(:), Z(:)] - dim) * dx / 2;
    
    % get fields
    data_raw{i} = fscanf(fs, '%e');
    data_raw{i} = reshape(data_raw{i}, 6, [], step);
    E{i} = data_raw{i}(1:3, :, 1:end-1);
    H{i} = data_raw{i}(4:6, :, :);
    H{i} = (H{i}(:, :, 1:end - 1) + H{i}(:, :, 2:end)) / 2;
    
    % process fourier transform
    E{i} = fft(E{i}, step - 1, 3) ./ reshape(ref_f, 1, 1, step - 1);
    H{i} = fft(H{i}, step - 1, 3) ./ reshape(ref_f, 1, 1, step - 1);
    
    fclose(fs);
end






