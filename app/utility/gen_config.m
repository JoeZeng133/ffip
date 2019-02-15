function config_str = gen_config(basic, medium, geometry, source, varargin)
%GEN_CONFIG Summary of this function goes here
%   Detailed explanation goes here
config_str = "";

% basic configuration
config_str = config_str + sprintf('basic {\n');
config_str = config_str + sprintf('%e %e\n', basic.dt, basic.dx);
config_str = config_str + sprintf('%d %d %d\n', basic.dim);
config_str = config_str + sprintf('%d\n', basic.sf_layer);
config_str = config_str + sprintf('%d\n', basic.tf_layer);
config_str = config_str + sprintf('%d\n', basic.PML_d);
config_str = config_str + sprintf('%d\n', basic.step);
config_str = config_str + sprintf('%e %e\n', basic.er_bg, basic.ur_bg);
config_str = config_str + sprintf('}\n');

% medium configuration
config_str = config_str + sprintf('medium %d {\n', numel(medium));
for i = 1 : numel(medium)
    item = medium{i};
    config_str = config_str + "{";
    config_str = config_str + sprintf('%e %e %e %e %d\n', item.er, item.sigma_e, item.ur, item.sigma_u, numel(item.poles));
    for j = 1 : numel(item.poles)
        pole = item.poles{j};
        pole_str = "";
        if pole.type == "Lorentz"
            pole_str = sprintf('{ Lorentz %e %e %e }\n', pole.rel_perm, pole.freq, pole.damp);
        end
        
        if pole.type == "Drude"
            pole_str = sprintf('{ Drude %e %e }\n', pole.freq, pole.inv_relaxation);
        end
        
        if pole.type == "Deybe"
            pole_str = sprintf('{ Deybe %e %e }\n', pole.rel_perm, pole.relaxation);
        end
        
        if pole.type == "CP"
            pole_str = sprintf('{ CP %e %e %e %e }\n',pole.A, pole.phi, pole.Omega, pole.Gamma);
        end
        
        config_str = config_str + pole_str;
    end
    config_str = config_str + sprintf('}\n');
end
config_str = config_str + sprintf('}\n');



% geometry
config_str = config_str + sprintf('geometry %d {\n', numel(geometry));
for i = 1 : numel(geometry)
    item = geometry{i};
    geo_str = "";
    if item.type == "sphere"
        geo_str = sprintf('{ sphere %d %e %e %e %e }\n', item.medium_idx, item.radius, item.position);
    end
    
    if item.type == "box"
        geo_str = sprintf('{ box %d %e %e %e %e %e %e }\n', item.medium_idx, item.lower_position, item.upper_position);
    end
    
    if item.type == "disk"
        geo_str = srpintf('{ disk %d %e %e %e %e %e %d }\n', item.medium_idx, item.position, item.radius, item.height, item.direction);
    end
    
    if item.type == "inhom"
        tmp = size(item.rho);
        if norm(tmp(:) - item.dim(:)) ~= 0
            error('Invalid dimension');
        end
        
        geo_str = sprintf('{ inhom %d %d %s }\n', item.medium_idx1, item.medium_idx2, item.filename);
        fileID = fopen(item.filename, 'w');
        fprintf(fileID, '%d %d %d\n', item.dim);
        for j = 1 : 3
            fprintf(fileID, '%e %e\n', item.position(j), item.dx(j));
        end
        fprintf(fileID, '%e ', item.rho(:));
        fclose(fileID);
    end
    
    if item.type == "blob"
        blob = item;
        config_str = config_str + sprintf('{ blob %s }\n', blob.input_file);
        fileID = fopen(blob.input_file, 'w');
        fprintf(fileID, '%d\n', numel(blob.x));
        fprintf(fileID, '%e %e %e %d %e\n', [blob.x(:), blob.y(:), blob.z(:), blob.medium_idx(:), blob.amp(:)]');
        fclose(fileID);
    end
    
    config_str = config_str + geo_str;
end
config_str = config_str + sprintf('}\n');

% source
config_str = config_str + sprintf('source %d {\n', numel(source));
for i = 1 : numel(source)
    item = source{i};
    s_str = "";
    if item.type == "plane"
        s_str = sprintf('{ plane %d %d %s %e %e %e }\n', item.dim_neg, item.dim_pos, item.func_type, item.fp, item.delay, item.ref_pos);
    end
    
    if item.type == "dipole"
        s_str = sprintf('{ dipole %s }\n', item.filename);
        output = [item.x(:), item.y(:), item.z(:), item.amp(:), item.func_type(:), item.fp(:), item.delay(:), item.ctype(:)];
        fileID = fopen(item.filename, 'w');
        fprintf(fileID, '%d\n', numel(item.x));
        fprintf(fileID, '%e %e %e %e %c %e %e %d\n', output');
        fclose(fileID);
    end
    
    config_str = config_str + s_str;
end
config_str = config_str + sprintf('}\n');

%optional configurations
p = inputParser();
addParameter(p, 'nearfield', 0);
addParameter(p, 'farfield', 0);
addParameter(p, 'flux', 0);
addParameter(p, 'nearfield_convergence', 0);
addParameter(p, 'nearfield_time', 0);
addParameter(p, 'config_filename', 'config.in');
addParameter(p, 'step_output', 0);
addParameter(p, 'num_proc', 0);

parse(p, varargin{:});
res = p.Results;

if isstruct(res.nearfield)
    nf = res.nearfield;
    config_str = config_str + sprintf('nearfield %s %s \n', nf.input_file, nf.output_file);
    fileID = fopen(nf.input_file, 'w');
    fprintf(fileID, '%d\n', numel(nf.x));
    fprintf(fileID, '%e %e %e %e\n', [nf.x(:), nf.y(:), nf.z(:), nf.freq(:)]');
    fclose(fileID);
end

if isstruct(res.farfield)
    ff = res.farfield;
    config_str = config_str + sprintf('farfield %s %s \n', ff.input_file, ff.output_file);
    fileID = fopen(ff.input_file, 'w');
    fprintf(fileID, '%e %e %e\n', ff.lower_position);
    fprintf(fileID, '%e %e %e\n', ff.upper_position);
    fprintf(fileID, '%d\n', numel(ff.theta));
    fprintf(fileID, '%e %e %e %e\n', [ff.theta(:), ff.phi(:), ff.rho(:), ff.freq(:)]');
    fclose(fileID);
end

if isstruct(res.flux)
    flux = res.flux;
    config_str = config_str + sprintf('flux %s %s \n', flux.input_file, flux.output_file);
    fileID = fopen(flux.input_file, 'w');
    fprintf(fileID, '%e %e %e\n', flux.lower_position);
    fprintf(fileID, '%e %e %e\n', flux.upper_position);
    fprintf(fileID, '%d\n', numel(flux.freq));
    fprintf(fileID, '%e\n', flux.freq(:));
    fclose(fileID);
end

if res.step_output == 0
    config_str = config_str + sprintf('Stop_Step_Output\n');
end

if res.num_proc ~= 0
    config_str = config_str + sprintf('num_proc %d\n', res.num_proc);
end

fileID = fopen(res.config_filename, 'w');
fprintf(fileID, '%s', config_str);
fclose(fileID);
end

