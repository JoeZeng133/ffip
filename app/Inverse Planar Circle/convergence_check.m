close
clc
clear

fileID = fopen('debug.out', 'r');
n = fscanf(fileID, '%d', 1);
data = fscanf(fileID, '%e');
fclose(fileID);
data = reshape(data, n * 6, []);
%%

for i = 1 : 20 : size(data, 1)
    plot(data(i, :)), hold on
end

hold off
