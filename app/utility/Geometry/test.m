s1 = make_disk([0 0], 2, 'num_elements', 20);
r1 = make_rect([0 0], [1 1]);
re = subtract(s1 ,r1);

[tfin, x, y, z, p1, p2] = make_inhom(re, 1, [-4 -4], [4 4], 0.1);


plot(re), hold on
plot(x(tfin), y(tfin), '.')
plot(x(~tfin), y(~tfin), 'x'), hold off
