%% manually read fields at each time and convert to N, L to far fields
clc
clear
close all

e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
c0 = 3e8;
eta0 = sqrt(u0 / e0);
er =  1;
e = e0 * er;
ur = 1;
u = u0;
c = 1 ./ sqrt(e * u);
eta = sqrt(u ./ e);


Sc = 1 / sqrt(3);
dt = 7.2e-17;
dx = c0 * dt / Sc;
dim = [10, 10, 10];
step = 600;

Np = 20;                            %center frequency of the rickerwavelet
fp = c0 / (Np * dx);
ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);

t = (1:step) * dt;
df = 1 / (step - 1) / dt;
ft = linspace(0.5 * fp, 1.5 * fp, 20)';
ref_signal = ricker(t, fp, 1/fp);
ref = sum(ref_signal .* exp(-1j * 2 * pi * ft * t), 2);
slice_len = numel(ft);
k = ft * 2 * pi / c0;

rho = linspace(1000, 2000, 11) * (Np * dx);
phi = linspace(0, 2 * pi, 10);
th = linspace(0 * pi,  pi, 10);

[Th, Phi, Rho] = ndgrid(th, phi, rho);
Th = Th(:);
Phi = Phi(:);
Rho = Rho(:);
[Xo, Yo, Zo] = sph2cart(Phi, pi / 2 - Th, 1);     %r0 observation directions

pos = [ft(:), ones([numel(ft), 3])];
fileID = fopen('request.in', 'w');
fprintf(fileID, '%d\n', size(pos, 1));
fprintf(fileID, '%e %e %e %e\n', pos');
fclose(fileID);

%% theoretical fields
[Th2, Phi2, Rho2, Ft2] = ndgrid(th, phi, rho, ft);
Th2 = Th2(:);
Phi2 = Phi2(:);
Rho2 = Rho2(:);
Ft2 = Ft2(:);
Omega = 2 * pi * Ft2;

K2 = Omega ./ c;

Hphi_p = 1 / (4 * pi) * exp(-1j * K2 .* Rho2) .* (1j * K2 ./ Rho2 + 1 ./ Rho2.^2) .* sin(Th2);
Er_p = 1 / (4 * pi) * exp(-1j * K2 .* Rho2) .* (2 * eta ./ Rho2.^2 + 2 ./ (1j * Omega .* e .* Rho2.^3)) .* cos(Th2);
Eth_p = 1 / (4 * pi) * exp(-1j * K2 .* Rho2) .* (1j * Omega .* u ./ Rho2 + 1 ./ (1j * Omega .* e .* Rho2.^3) + eta ./ Rho2.^2) .* sin(Th2);


%% Numrical Fields
% Input Format: each N2F face has an output file "face[i].out"
% in each "face[i].out"
% 1: p1 (lower corner coordinate in computational frame)
% 2: p2 (upper corner coordinate in computational frame)
% 3: normal vector
% 5...: Jx Jy Jz Mx My Mz at each (position, frequency). 
% Positions are ordered in the
% way that abides to the local reference frame of each face (x1, x2, x3).
% where x3 is the normal direction In this way, x1 is incremented first,
% then x2, then x3
clc

data_raw = cell([6 1]);
int_weights = cell([6 1]);      %N2F integration weights
pos = cell([6 1]);              %N2F sample positions
make_complex = @(x, y) x + 1j * y;

Ntot = 0;
Ltot = 0;
for i = 1 : 6
    
    filename = ['face', num2str(i - 1), '.out'];
    fs = fopen(filename);
    p1 = fscanf(fs, '%e %e %e', [3 1]);
    p2 = fscanf(fs, '%e %e %e', [3 1]);
    norm_vec = fscanf(fs, '%e %e %e', [3 1]);
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
    
    if (norm_vec(1) ~= 0)
        [inty, intz] = ndgrid(int{2}, int{3});
        int_weights{i} = inty .* intz;
        int_weights{i} = int_weights{i}(:) * dx^2;
        [Y, Z, X] = ndgrid(p1(2):2:p2(2), p1(3):2:p2(3), p1(1):2:p2(1));
    else
        if (norm_vec(2) ~= 0)
            [intz, intx] = ndgrid(int{3}, int{1});
            int_weights{i} = intz .* intx ;
            int_weights{i} = int_weights{i}(:) * dx^2;
            [Z, X, Y] = ndgrid(p1(3):2:p2(3), p1(1):2:p2(1), p1(2):2:p2(2));
        else
            [intx, inty] = ndgrid(int{1}, int{2});
            int_weights{i} = intx .* inty ;
            int_weights{i} = int_weights{i}(:) * dx^2;
            [X, Y, Z] = ndgrid(p1(1):2:p2(1), p1(2):2:p2(2), p1(3):2:p2(3));
        end
    end
        
    
    
    % get sample position - center for each point
    pos{i} = ([X(:), Y(:), Z(:)] - dim) * dx / 2;
    
    % get fields
    data_raw{i} = fscanf(fs, '%e');
    data_raw{i} = reshape(data_raw{i}, 2, []);
    data_raw{i} = make_complex(data_raw{i}(1, :), data_raw{i}(2, :));
    data_raw{i} = reshape(data_raw{i}, 6, [], slice_len);
    J = data_raw{i}(1:3, :, :) ./ reshape(ref, 1, 1, []);             % 3 x pos x slice_len
    M = data_raw{i}(4:6, :, :) ./ reshape(ref, 1, 1, []);
  
    % get N, L
    arg = reshape([Xo, Yo, Zo], [], 1, 3) .* reshape(pos{i}, 1, [], 3); %ro * r', obs x pos
    arg = sum(arg, 3);  
    arg = exp(1j * arg .* reshape(k, 1, 1, [])); %exp(j * ro * r' * k), obs x pos x slice_len
    arg = arg .* reshape(int_weights{i}, 1, []); %exp(j * ro * r' * k)ds, obs x pos x slice_len
    
    % Jexp(j * ro * r' * k)ds, 3 x obs x pos x slice_len
    N = reshape(J, 3, 1, [], slice_len) .* reshape(arg, [1 size(arg)]);
    L = reshape(M, 3, 1, [], slice_len) .* reshape(arg, [1 size(arg)]);
    N = reshape(sum(N, 3), 3, [], slice_len);   %3 x obs x slice_len
    L = reshape(sum(L, 3), 3, [], slice_len);
    
    Ntot = Ntot + N;
    Ltot = Ltot + L;
    fclose(fs);
end

% get Nth, Lth From N, L
proj_th = [cos(Th) .* cos(Phi), cos(Th) .* sin(Phi), -sin(Th)]';   %3 x obs
proj_phi = [-sin(Phi), cos(Phi), zeros(size(Phi))]';                 
Nth = reshape(sum(proj_th .* Ntot, 1), [], slice_len);                  %obs x slice_len
Lth = reshape(sum(proj_th .* Ltot, 1), [], slice_len);
Nphi = reshape(sum(proj_phi .* Ntot, 1), [], slice_len);
Lphi = reshape(sum(proj_phi .* Ltot, 1), [], slice_len);

% get Eth, Hth, Ephi, Hphi
arg = (1j ./ (4 * pi * Rho) * k') .* (exp(-1j * Rho * k'));         %obs x slice_len
Eth = -(Lphi + eta0 * Nth) .* arg;
Ephi = (Lth - eta0* Nphi) .* arg;
Hth = (Nphi - Lth / eta0) .* arg;
Hphi = -(Nth + Lphi / eta0) .* arg; 

%% correlation plot
% it should be noted that the phase is not very accurant so is left out
% while the absolute value has a pleasant agreement with theoretical
% values.

figure(3)
plot(abs(Eth(:)), abs(Eth_p(:)), '.')
axis equal

figure(4)
plot(abs(Hphi(:)), abs(Hphi_p(:)), '.')
axis equal
