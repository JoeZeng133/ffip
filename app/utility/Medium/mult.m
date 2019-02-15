function [m] = mult(m, f)
%MULT Summary of this function goes here
%   Detailed explanation goes here
m.er = m.er * f;
m.sigma_e = m.sigma_e * f;
m.ur = m.ur * f;
m.sigma_u = m.sigma_u * f;

for i = 1 : numel(m.poles)
    pole = m.poles{i};
    if strcmp(pole.type, 'Drude')
        m.poles{i}.freq = m.poles{i}.freq * sqrt(f);
    end
    
    if strcmp(pole.type, 'Lorentz')
        m.poles{i}.rel_perm = m.poles{i}.rel_perm * f;
    end
    
    if strcmp(pole.type, 'CP')
        m.poles{i}.A = m.poles{i}.A * f;
    end
end
end

