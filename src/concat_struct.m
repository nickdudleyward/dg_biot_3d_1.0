function out = concat_struct(varargin)

% Concatenate struct vectors, allowing different fields in each vector.
% Arguments should be 1-dimensional struct arrays, rows or columns.
% Output is a 1-dimensional struct array, either a row or a column
% according to the first argument (a row if the first argument is a scalar, 
% an empty struct if there are no arguments)
% Concatenation is performed along the length of each input array, adding 
% new fields to the new struct array as required.
% SPE 2019/12/01

if nargin == 0
    out = struct();
    return
end

out = varargin{1};
k   = length(out);

for j = 2:nargin
    s = varargin{j};
    if ~isstruct(s) || ~isvector(s)
        error('All arguments must be 1-dimensional struct arrays')
    end
    f = fieldnames(s);
    for n=1:length(s)
        k = k+1;
        for fn = 1:length(f)
            out(k).(f{fn}) = s(n).(f{fn});
        end
    end
end



