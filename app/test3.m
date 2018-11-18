% geometry test
clear
clc

data = load('data.in');
x = -1:0.01:1-0.01;
y = -1:0.01:1-0.01;

data = reshape(data, 200, 200);
imagesc(x, y, data)
xlabel('x')
ylabel('y')