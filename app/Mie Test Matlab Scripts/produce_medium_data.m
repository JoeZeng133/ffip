clear
clc

freq_span = linspace(300e12, 1500e12, 100)';
lam = 3e8 ./ freq_span;
er = conj(Au(2 * pi * freq_span));
rel_er = real(er);
loss_tangent = imag(er) ./ rel_er;

plot(lam, real(er)), hold on
plot(lam, imag(er)), hold off

