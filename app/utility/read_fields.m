function varargout = read_fields(filename, varargin)
%READ_FIELDS [E1, H1, E1, H2, ...] = read_fields(filename, a1, a2, ...)
%   Detailed explanation goes here
data = load(filename);
make_complex = @(x, y) x + 1j * y;

p = inputParser;
addParameter(p, 'partition', 0);
addParameter(p, 'norm', 1);
parse(p, varargin{:});

E = [make_complex(data(:, 1), data(:, 2)), make_complex(data(:, 3), data(:, 4)), make_complex(data(:, 5), data(:, 6))];
H = [make_complex(data(:, 7), data(:, 8)), make_complex(data(:, 9), data(:, 10)), make_complex(data(:, 11), data(:, 12))];
E = E ./ p.Results.norm(:);
H = H ./ p.Results.norm(:);

if numel(p.Results.partition) == 1
    varargout{1} = E;
    varargout{2} = H;
else
    start = 1;
    partition = p.Results.partition;
    if sum(partition) ~= size(E, 1)
        error('Invalid Partitions')
    end
    
    for i = 1 : numel(partition)
        varargout{2 * i - 1}    = E(start:start + partition(i) - 1, :); 
        varargout{2 * i}        = H(start:start + partition(i) - 1, :);
        start = start + partition(i);
    end
end
end

