val_real = h5read('test_output.h5', '/val_real');
val_imag = h5read('test_output.h5', '/val_imag');
val = val_real + 1j * val_imag;

%% 

for i = 1 : 5 : size(val, 1)
    figure(i)
    surf(abs(squeeze(val(i, :, :))))
    title(num2str(i))
end