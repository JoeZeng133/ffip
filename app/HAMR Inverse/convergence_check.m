close
clc
clear

file = fopen('slice.out', 'r');
p1 = fscanf(file, '%e %e %e', [1 3]);
p2 = fscanf(file, '%e %e %e', [1 3]);
dim = (p2 - p1) / 2 + 1;
time_step = 3000;
data = fscanf(file, '%e');
fclose(file);

%%
data = reshape(data, [dim(2) dim(3) time_step]);
norm = max(abs(data(:)));
timeinterval = 0.005;
n = size(data, 1);
figure
tic
for i = 1 : 10 : time_step
    st = toc;
    surf(data(:, :, i) / norm)
    axis([-inf inf -inf inf -1 1])
    xlabel('x')
    ylabel('y')
    while(toc - st < timeinterval)
    end
    getframe;
end