x = linspace(0, 1, 5);
y = linspace(0, 1, 10);
[X, Y] = ndgrid(x, y);
data = X.^2 + Y;

h5create('test.h5', '/source/sth', size(data), 'Datatype', 'double');
h5write('test.h5', '/source/sth', data);
h5writeatt('test.h5', '/source/sth', 'x', x);
h5writeatt('test.h5', '/source/sth', 'y', y);
