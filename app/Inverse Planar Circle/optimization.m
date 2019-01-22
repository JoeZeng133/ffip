clear
clc
% close all

load('case_configuration.mat');
rho = ones(inhom_dim(1:2) + 1);                          %start with rho=1
tol = 1e-6;
step = 0.1; 
a = 0.9;
max_itr = 30;
itr = 0;
max_size = [1 1000];

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
        surf(rho)
%         pcolor(rho)
        colorbar
        shading flat
        caxis([0 1])
        zlim([0 1])
        title(['Iteration ', num2str(itr)])
        xlabel('y')
        ylabel('x')
       
        getframe;
    end
    
    rho_list{itr} = rho;
    % forward simulation configuration
    % write to geometry file
    rho_out = rho .* ones([1, 1, inhom_dim(3) + 1]);
    file_geometry_forward = 'forward_geometry.in';
    fileID = fopen(file_geometry_forward, 'w');
    fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
    for i = 1 : 3
        fprintf(fileID, '%e %e\n', inhom_p1(i), inhom_dx(i));
    end
    fprintf(fileID, '%e\n', rho_out(:));
    fclose(fileID);

    % write to probes file
    file_probes_input_forward = 'forward_probes.in';
    fileID = fopen(file_probes_input_forward, "w");
    fprintf(fileID, "%d\n", num_probes + num_inhom_pos);
    fprintf(fileID, "%e %e %e %e\n", [[probes_pos; inhom_pos], ft * ones([num_probes + num_inhom_pos, 1])]');
    fclose(fileID);

    % basic configuration
    fileID = fopen('config.in', 'w');
    fprintf(fileID, 'basic {\n');
    fprintf(fileID, '%e %e\n', dt, dx);
    fprintf(fileID, '%d %d %d\n', dim);
    fprintf(fileID, '%d\n', sf_layer);
    fprintf(fileID, '%d\n', tf_layer);
    fprintf(fileID, '%d\n', PML_d);
    fprintf(fileID, '%d\n', time_step);
    fprintf(fileID, '%e %e\n', er_bg, ur_bg);
    fprintf(fileID, '}\n');

    % medium configuration
    fprintf(fileID, "medium 3 {\n");
    % medium 0, background medium
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
    fprintf(fileID, " }\n");
    % medium 1 Au
    Au_print(fileID);
    % medium 2 FePt
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", real(fept_er), fept_sig, ur_bg, 0);
    fprintf(fileID, " }\n");

    fprintf(fileID, "}\n");

    % geometry configuration
    fprintf(fileID, "geometry 3 {\n");
    % geometry 0 gold aperture
    fprintf(fileID, "{ ");
    fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_forward);
    fprintf(fileID, " }\n");
    % geometry 1 Gold layer
    fprintf(fileID, "{ ");
    fprintf(fileID, "box %d %e %e %e %e %e %e", 1, [0, 0, inhom_p1(3)], [p2(1), p2(2), inhom_p2(3)]);
    fprintf(fileID, " }\n");
    % % geometry 2 FePt Layer
    fprintf(fileID, "{ ");
    fprintf(fileID, "box %d %e %e %e %e %e %e", 2, fept_p1, fept_p2);
    fprintf(fileID, " }\n");

    fprintf(fileID, "}\n");

    % plane wave source
    fprintf(fileID, "source 1 {\n");
    fprintf(fileID, "{ ");
    fprintf(fileID, 'plane %d %d %s %e %e %e', projector_padding, dim(3) + projector_padding, 's', fs, delay, 0);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");
    
    % no step number output
    fprintf(fileID, 'Stop_Step_Output\n');
    
    % probes
    file_probes_output_forward = 'forward_output.out';
    fprintf(fileID, "nearfield %s %s\n", file_probes_input_forward, file_probes_output_forward);
    fclose(fileID);

    save("forward_results");
    disp('forward simulation configuration complete');
    % forward simulation results
    call_exe('std_config')
    ref_plane = load('reference.out');
    ref_plane = reshape(ref_plane, 1, []);
    t = (0:numel(ref_plane) - 1) * dt;
    ref_plane_fft = sum(ref_plane .* exp(-1j * 2 * pi * ft * t), 2);

    data = load(file_probes_output_forward);

    E_forward = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
    H_forward = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
    E_forward = E_forward / ref_plane_fft;
    H_forward = H_forward / ref_plane_fft;

    E_forward_probes = E_forward(1:num_probes, :);
    H_forward_probes = H_forward(1:num_probes, :);
    E_forward_inhom = E_forward(num_probes + 1:end, :);
    H_forward_inhom = H_forward(num_probes + 1:end, :);

%     target_func = 0.5 * sum(abs(E_forward_probes(:, 3) - E_target_probes(:, 3)).^2, 1); %only Ez
    obj_func(itr) = 0.5 * sum(sum(abs(H_forward_probes - H_target_probes).^2, 2), 1);
    fprintf('Current objective function is evaluated at : %e\n', obj_func(itr));
    
    % adjoint simulation configuration
    ricker = @(t, fp, d) (1 - 2 * (pi * fp * (t - d)).^2) .* exp(-(pi * fp * (t - d)).^2);
    ref_dipole = sin(2 * pi * fs * t);
    ref_dipole_fft = sum(ref_dipole.* exp(-1j * 2 * pi * ft * t), 2);
    
    Ex = 1;
    Ey = 2;
    Ez = 4;
    Hx = 7 - Ex;
    Hy = 7 - Ey;
    Hz = 7 - Ez;
    dipoles = -conj(H_forward_probes - H_target_probes);
    % dipoles = 2 * (abs(E_forward_probes) - abs(E_target_probes)) .* (E_forward_probes) ./ abs(E_forward_probes);
    amp = abs(dipoles);
    dipoles_delay =  - (angle(dipoles))  / (2 * pi * ft);
    ctype = [Hx * ones([num_probes, 1]); Hy * ones([num_probes, 1]); Hz * ones([num_probes, 1])];
    dipoles_pos = [probes_pos;probes_pos;probes_pos];
    
%     dipoles = conj(E_forward_probes(:, 3) - E_target_probes(:, 3));
%     amp = abs(dipoles);
%     dipoles_delay = delay - angle(dipoles)  / (2 * pi * ft);
%     ctype = Ez * ones([num_probes, 1]);
%     dipoles_pos = probes_pos * dx;
    
    % write dipoles
    
    output_dipoles = [dipoles_pos, amp(:), 's' * ones(size(amp(:))), fs * ones(size(amp(:))), dipoles_delay(:), ctype];
    file_dipoles_adjoint = 'adjoint_dipoles.in';
    fileID = fopen(file_dipoles_adjoint, 'w');
    fprintf(fileID, '%d\n', size(output_dipoles, 1));
    fprintf(fileID, '%e %e %e %e %s %e %e %d\n', output_dipoles');
    fclose(fileID);

    % basic configuration
    fileID = fopen('config.in', 'w');
    fprintf(fileID, 'basic {\n');
    fprintf(fileID, '%e %e\n', dt, dx);
    fprintf(fileID, '%d %d %d\n', dim);
    fprintf(fileID, '%d\n', sf_layer);
    fprintf(fileID, '%d\n', tf_layer);
    fprintf(fileID, '%d\n', PML_d);
    fprintf(fileID, '%d\n', time_step);
    fprintf(fileID, '%e %e\n', er_bg, ur_bg);
    fprintf(fileID, '}\n');

    % medium configuration
    fprintf(fileID, "medium 3 {\n");
    % medium 0, background medium
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
    fprintf(fileID, " }\n");
    % medium 1 Au
    Au_print(fileID);
    % medium 2 FePt
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", real(fept_er), fept_sig, ur_bg, 0);
    fprintf(fileID, " }\n");

    fprintf(fileID, "}\n");

    % geometry configuration
    fprintf(fileID, "geometry 3 {\n");
    % geometry 0 gold aperture
    fprintf(fileID, "{ ");
    fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_forward);
    fprintf(fileID, " }\n");
    % geometry 1 Gold layer
    fprintf(fileID, "{ ");
    fprintf(fileID, "box %d %e %e %e %e %e %e", 1, [0, 0, inhom_p1(3)], [p2(1), p2(2), inhom_p2(3)]);
    fprintf(fileID, " }\n");
    % % geometry 2 FePt Layer
    fprintf(fileID, "{ ");
    fprintf(fileID, "box %d %e %e %e %e %e %e", 2, fept_p1, fept_p2);
    fprintf(fileID, " }\n");

    fprintf(fileID, "}\n");

    % dipole sources
    fprintf(fileID, "source 1 {\n");
    fprintf(fileID, "{ ");
    fprintf(fileID, "dipole %s", file_dipoles_adjoint);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");
    
    % no step number output
    fprintf(fileID, 'Stop_Step_Output\n');
    
    % probes
    file_probes_output_adjoint = 'adjoint_output.out';
    fprintf(fileID, "nearfield %s %s\n", file_probes_input_forward, file_probes_output_adjoint );
    fclose(fileID);

    disp('adjoint simulation configuration complete');
    % adjoint fields
    call_exe('std_config')
    data = load(file_probes_output_adjoint);
    E_adjoint = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
    H_adjoint = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
    E_adjoint = E_adjoint / ref_dipole_fft;
    H_adjoint = H_adjoint / ref_dipole_fft;

    E_adjoint_probes = E_adjoint(1:num_probes, :);
    H_adjoint_probes = H_adjoint(1:num_probes, :);
    E_adjoint_inhom = E_adjoint(num_probes + 1:end, :);
    H_adjoint_inhom = H_adjoint(num_probes + 1:end, :);

    save("adjoint_results");
    disp('adjoint simulation results saved');
    % Sensitivity Calculation
    int_weights = @(dim) [0.5, ones([1, dim - 2]), 0.5];
    [w_X, w_Y, w_Z] = ndgrid(int_weights(inhom_dim(1) + 1), int_weights(inhom_dim(2) + 1), int_weights(inhom_dim(3) + 1));
    w = w_X .* w_Y .* w_Z;
    w = w(:) * dx^3;
    Se = 1j * (2 * pi * ft) * sum(E_forward_inhom .* E_adjoint_inhom, 2);
    de_drho = e0 * (er_const - er_bg);
    A = real(Se .* w * de_drho);
    A = reshape(A, inhom_dim + 1);
    A = sum(A, 3);

    Se_list{itr} = Se;
    sens{itr} = A;
    
    disp('sensitivity calculation complete')
    % Optimization
    step = a * step;
    
    lb = max(-rho(:), -step);
    ub = min(1 - rho(:), step);
    delta_rho = linprog(A(:), [], [], [], [], lb, ub);
    
%     delta_rho = -A(:) * step;
    
    exp_obj_func(itr) = obj_func(itr) + sum(delta_rho .* A(:));
    fprintf('The new objective function is expected to be : %e\n', exp_obj_func(itr));

    % updates
    rho = reshape(delta_rho, size(rho)) + rho;
    disp('updated');
    
end




