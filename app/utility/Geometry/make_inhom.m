function [tfin, x, y, z, p1, p2, cell_len] = make_inhom(shape, height, bbox_p1, bbox_p2, dx)
p1 = [bbox_p1(:); -height/2];
p2 = [bbox_p2(:); +height/2];

dim = round(p2 - p1) / dx;
cell_len = (p2 - p1) ./ dim;
dim = dim + 1;
[x, y, z] = ndgrid(...
    linspace(p1(1), p2(1), dim(1)),...
    linspace(p1(2), p2(2), dim(2)),...
    linspace(p1(3), p2(3), dim(3)));

[tx, ty] = ndgrid(linspace(p1(1), p2(1), dim(1)), linspace(p1(2), p2(2), dim(2)));
tfin = isinterior(shape, tx(:), ty(:));
tfin = reshape(tfin, size(tx));
tfin = repmat(tfin, 1, 1, dim(3));
end

