clear
clc
close all

data = load('projector.out');

%%

norm = max(abs(data(:)));
timeinterval = 0.005;
n = size(data, 1);
tic
for i = 1 : 1 : n
    st = toc;
    ex = data(i, :) / norm;
    plot(ex)
    axis([-inf inf -1 1])
    
    while(toc - st < timeinterval)
    end
    getframe;
end