%% testing utility->integral_ndim and interp_ndim
clear
clc

n = 10;
m = 1;
p = 5;
dx = 0.3;
dy = 0.11;
dz = 0.34;

x = (0:n-1) * dx;
y = (0:m-1) * dy;
z = (0:p-1) * dz;

[X, Y, Z] = ndgrid(x, y, z);

V = sin(X) + X.^2 + sin(Y) + Y.^2; + sin(Z) + Z.^2;

f = fopen('data.in', 'w');

fprintf(f, '%d %d %d\n', n, m, p);
fprintf(f, '%f %f\n', x(1), dx);
fprintf(f, '%f %f\n', y(1), dy);
fprintf(f, '%f %f\n', z(1), dz);
fprintf(f, '%f ', V(:));



fclose(f);
%%
f2 = fopen('request.in', 'w');
r = 10;
Xq = rand([1 r]) * x(end);
Yq = rand([1 r]) * y(end);
Zq = rand([1 r]) * z(end);

fprintf(f2, '\n%d\n', r);
fprintf(f2, '%f %f %f\n', [Xq;Yq;Zq]);
fclose(f2);



