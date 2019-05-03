m_real = h5read('mie_material.h5', '/m_real');
m_imag = h5read('mie_material.h5', '/m_imag');
x = h5read('mie_material.h5', '/x');
l = h5read('mie_material.h5', '/l');
m = m_real - 1j * m_imag;
num = numel(x);

qsca = zeros([num 1]);
qabs = zeros([num 1]);

for i = 1 : num
    res = Mie(m(i), x(i));
    qsca(i) = res(2);
    qabs(i) = res(3);
end

subplot(1, 2, 1)
plot(l, qsca)
subplot(1, 2, 2)
plot(l, qabs)