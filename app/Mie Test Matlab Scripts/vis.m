clear
clc
close all

file = fopen('debug.out', 'r');
dim = fscanf(file, '%d %d %d', [1 3]);
data = fscanf(file, '%e');
data = reshape(data, dim(1), dim(2), dim(3), []);
fclose(file);

%%

norm = max(abs(data(:)));
timeinterval = 0.005;
n = size(data, 4);
tic
for i = 1 : 1 : n
    st = toc;
    ex = reshape(data(:, :, :, i), dim(2:3)) / norm;
    surf(ex)
    axis([-inf inf -inf inf -1 1])
    
    while(toc - st < timeinterval)
    end
    getframe;
end