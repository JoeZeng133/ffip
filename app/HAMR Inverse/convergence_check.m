close
clc
clear

data = load('debug.out');

%%

for i = 1 : 10 : size(data, 2)
    plot(data(:, i)), hold on
end

hold off
