function [tfin, x, y, z, p1, p2, dim, dx] = make_inhom(shape, height, bbox_p1, bbox_p2, varargin)
% shape, height, p1, p2, dx, float(true, false), center([0 0 0])
p = inputParser;
addParameter(p, 'center', [0 0 0]);
addParameter(p, 'float', false);
addParameter(p, 'dx', [0 0 0]);
addParameter(p, 'dim', [0 0 0]);
parse(p, varargin{:});

p1 = [bbox_p1(:); -height/2];
p2 = [bbox_p2(:); +height/2];

[x, y, z, dim, dx] = make_interp_region(p1, p2, 'dim', p.Results.dim, 'dx', p.Results.dx);

[tx, ty] = ndgrid(linspace(p1(1), p2(1), dim(1)), linspace(p1(2), p2(2), dim(2)));
tfin = isinterior(shape, tx(:), ty(:));
tfin = reshape(tfin, size(tx));
tfin = repmat(tfin, 1, 1, dim(3));

center = p.Results.center;
x = x + center(1);
y = y + center(2);
z = z + center(3);
p1 = p1 + center(:);
p2 = p2 + center(:);

if p.Results.float
    tfin = tfin * 1;
end
end

