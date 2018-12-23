clear
clc
close all

load('case_configuration.mat');
rho = ones(inhom_dim(1:2) + 1);                          %start with rho=1
tol = 1e-4;
step = 0.2;
a = 0.9;
max_itr = 100;
itr = 0;

%% optimization
while (itr < max_itr)
    itr = itr + 1;
    fprintf('################# Iteration %d\n', itr);
    % forward simulation configuration
    % write to geometry file
    rho_out = rho .* ones([1, 1, inhom_dim(3) + 1]);
    file_geometry_forward = 'forward_geometry.in';
    fileID = fopen(file_geometry_forward, 'w');
    fprintf(fileID, '%d %d %d\n', inhom_dim + 1);
    fprintf(fileID, '%e %e\n', inhom_p1(1) * dx, dx);
    fprintf(fileID, '%e %e\n', inhom_p1(2) * dx, dx);
    fprintf(fileID, '%e %e\n', inhom_p1(3) * dx, dx);
    fprintf(fileID, '%e ', rho_out);
    fclose(fileID);

    % write to probes file
    file_probes_input_forward = 'forward_probes.in';
    fileID = fopen(file_probes_input_forward, "w");
    fprintf(fileID, "%d\n", num_probes + num_inhom_pos);
    fprintf(fileID, "%e %e %e %e\n", [[probes_pos; inhom_pos] * dx, fp * ones([num_probes + num_inhom_pos, 1])]');
    fclose(fileID);

    % basic configuration
    fileID = fopen('config.in', 'w');
    fprintf(fileID, "basic {\n");
    fprintf(fileID, "%e %e\n", dt, dx);
    fprintf(fileID, "%d %d %d\n", dim);
    fprintf(fileID, "%d\n", time_step);
    fprintf(fileID, "%e %e\n", er_bg, ur_bg);
    fprintf(fileID, "%d\n", PMl_d);
    fprintf(fileID, "}\n");

    % medium configuration
    fprintf(fileID, "medium 2 {\n");
    % medium 0, background medium
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
    fprintf(fileID, " }\n");

    % medium 1, scatterer medium
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", real(er_const), sig_const, ur_bg, 0);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");

    % geometry configuration
    fprintf(fileID, "geometry 1 {\n");
    % geometry 0, the inhomogeneous region with mixed medium1 and medium0
    fprintf(fileID, "{ ");
    fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_forward);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");

    % plane wave source
    fprintf(fileID, "source 1 {\n");
    fprintf(fileID, "{ ");
    fprintf(fileID, "eigen %d %e %e", dim(3), fp, delay);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");

    % probes
    file_probes_output_forward = 'forward_output.out';
    fprintf(fileID, "probe %s %s\n", file_probes_input_forward, file_probes_output_forward);
    fclose(fileID);

    save("forward_results");
    disp('forward simulation configuration complete');
    % forward simulation results
    !./std_config
    data = load(file_probes_output_forward);

    E_forward = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
    H_forward = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
    E_forward = E_forward / ref_signal_fft;
    H_forward = H_forward / ref_signal_fft;

    E_forward_probes = E_forward(1:num_probes, :);
    H_forward_probes = H_forward(1:num_probes, :);
    E_forward_inhom = E_forward(num_probes + 1:end, :);
    H_forward_inhom = H_forward(num_probes + 1:end, :);

    target_func = 0.5 * sum(sum(abs(E_forward_probes - E_target_probes).^2, 2), 1);
    fprintf('Current objective function is evaluated at : %e\n', target_func);
    if (target_func < tol)
        break;
    end
    
    % adjoint simulation configuration
    Ex = 1;
    Ey = 2;
    Ez = 4;
    dipoles = conj(E_forward_probes - E_target_probes);
    % dipoles = 2 * (abs(E_forward_probes) - abs(E_target_probes)) .* (E_forward_probes) ./ abs(E_forward_probes);
    amp = abs(dipoles);
    dipoles_delay = delay - angle(dipoles)  / (2 * pi * fp);
    ctype = [Ex * ones([num_probes, 1]), Ey * ones([num_probes, 1]), Ez * ones([num_probes, 1])];
    dipoles_slice = (amp > 1e-3);
    dipoles_pos = [probes_pos;probes_pos;probes_pos] * dx;
    
    % write dipoles
    amp = amp(dipoles_slice);
    dipoles_delay = dipoles_delay(dipoles_slice);
    ctype = ctype(dipoles_slice);
    dipoles_pos = dipoles_pos(dipoles_slice(:), :);
    output_dipoles = [dipoles_pos, amp, fp * ones(size(amp)), dipoles_delay, ctype];
    file_dipoles_adjoint = 'adjoint_dipoles.in';
    fileID = fopen(file_dipoles_adjoint, 'w');
    fprintf(fileID, '%d\n', size(output_dipoles, 1));
    fprintf(fileID, '%e %e %e %e %e %e %d\n', output_dipoles');
    fclose(fileID);

    % basic configuration
    fileID = fopen('config.in', 'w');
    fprintf(fileID, "basic {\n");
    fprintf(fileID, "%e %e\n", dt, dx);
    fprintf(fileID, "%d %d %d\n", dim);
    fprintf(fileID, "%d\n", time_step);
    fprintf(fileID, "%e %e\n", er_bg, ur_bg);
    fprintf(fileID, "%d\n", PMl_d);
    fprintf(fileID, "}\n");

    % medium configuration
    fprintf(fileID, "medium 2 {\n");
    % medium 0, background medium
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", er_bg, 0, ur_bg, 0);
    fprintf(fileID, " }\n");

    % medium 1, scatterer medium
    fprintf(fileID, "{ ");
    fprintf(fileID, "%e %e %e %e 0", real(er_const), sig_const, ur_bg, 0);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");

    % geometry configuration
    fprintf(fileID, "geometry 1 {\n");
    % geometry 0, the inhomogeneous region with mixed medium1 and medium0
    fprintf(fileID, "{ ");
    fprintf(fileID, "inhom %d %d %s", 1, 0, file_geometry_forward);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");

    % dipole sources
    fprintf(fileID, "source 1 {\n");
    fprintf(fileID, "{ ");
    fprintf(fileID, "dipole %s", file_dipoles_adjoint);
    fprintf(fileID, " }\n");
    fprintf(fileID, "}\n");

    % probes
    file_probes_output_adjoint = 'adjoint_output.out';
    fprintf(fileID, "probe %s %s\n", file_probes_input_forward, file_probes_output_adjoint );
    fclose(fileID);

    disp('adjoint simulation configuration complete');
    % adjoint fields
    !./std_config
    data = load(file_probes_output_adjoint);
    E_adjoint = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
    H_adjoint = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
    E_adjoint = E_adjoint / ref_signal_fft;
    H_adjoint = H_adjoint / ref_signal_fft;

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
    Se = 1j * (2 * pi * fp) * sum(E_forward_inhom .* E_adjoint_inhom, 2);
    de_drho = e0 * (er_const - er_bg);
    A = real(Se .* w * de_drho);
    A = reshape(A, inhom_dim + 1);
    A = sum(A, 3);

    disp('sensitivity calculation complete')
    % Optimization
    step = a * step;
    lb = max(-rho(:), -step);
    ub = min(1 - rho(:), step);
    delta_rho = linprog(A(:), [], [], [], [], lb, ub);
    new_target_func = target_func + sum(delta_rho .* A(:));
    fprintf('The new objective function is expected to be : %e\n', new_target_func);

    % updates
    rho = reshape(delta_rho, size(rho)) + rho;
    disp('updated');
    if (mod(itr - 1, 3) == 0)
        figure
        pcolor(rho)
        colorbar
        shading flat
        title(['Iteration ', num2str(itr)])
        xlabel('y')
        ylabel('x')
        getframe;
    end
end