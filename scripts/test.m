[x, y, z] = read_geom('forwad_output.h5', '/geometry info request 0 Ey');


plot3(x, y, z, '.')
xlabel('x')
ylabel('y')
zlabel('z')