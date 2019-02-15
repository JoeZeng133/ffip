function m = add(m1, m2)
m.er = m1.er + m2.er;
m.sigma_e = m1.sigma_e + m2.sigma_e;
m.ur = m1.ur + m2.ur;
m.sigma_u = m1.sigma_u + m2.sigma_u;
m.poles = {m1.poles{:}, m2.poles{:}};
end