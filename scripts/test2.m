[x, y, z, f, v] = read_dft('gbeam_output.h5', '/volume fields dft 1');

v = squeeze(v);

surf(abs(v))