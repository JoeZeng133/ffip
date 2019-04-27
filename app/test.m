data1 = h5read('output2.h5', '/volume fields dft0/real');
data2 = h5read('output2.h5', '/volume fields dft0/imag');
data = data1 + 1j * data2;



%%
figure
plot(data1(:)), hold on
plot(data2(:)), hold on