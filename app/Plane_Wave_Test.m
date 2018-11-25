% plane wave projector test
clear
clc
close all

d_size = [600, 25, 25];
data = load('data.out');
data = reshape(data, d_size);
norm = max(abs(data(:)));

timeinterval = 0.005;
n = size(data, 1);
tic
for i = 1 : n
    st = toc;
    ex = reshape(data(i, :, :), d_size(2:3)) / norm;
    surf(ex)
    axis([-inf inf -inf inf -1 1])
    
    while(toc - st < timeinterval)
    end
    getframe;
end
% plot(data(:), 'r')
% plot(data(:, 2), 'b')