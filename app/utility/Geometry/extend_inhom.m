function [inhom] = extend_inhom(inhom, extension, padval)
extension = extension(:)';
if size(extension, 2) ~= 3
	error('extension size is incorrect');
end

inhom.rho = padarray(inhom.rho, extension, padval, 'both');
inhom.position = inhom.position(:) - inhom.dx(:) .* extension(:);
inhom.dim = size(inhom.rho);
inhom.p1 = inhom.position;
inhom.p2 = inhom.p1(:) + inhom.dx(:) .* (inhom.dim(:) - 1);
end