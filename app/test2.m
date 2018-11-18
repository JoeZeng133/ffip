% plane wave projector test

clear
clc

data = load('data.in');
n = size(data, 1);

timeinterval = 0.005;
tic
for i = 1 : n
    st = toc;
    ez = data(i, :);
    plot(ez)
    axis([-inf inf -1 1])
    
    while(toc - st < timeinterval)
    end
    getframe;
end