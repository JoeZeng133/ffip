r0 = 30;
a = 2;

f = @(x) -1/pi * atan(a * (x - r0)) + 0.5;
fd = @(x) a / pi ./ (1 + (a * (x - r0)).^2);

x = linspace(r0 - 10, r0 + 10, 100);

plot(x, f(x))
