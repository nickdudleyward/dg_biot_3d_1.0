function h5save(h5_file, datasetname, array, datatype)

% h5create, h5write insist on absolute names; add a leading / if not
% already present

if datasetname(1) ~= '/'
    datasetname = strcat('/', datasetname);
end

% if datatype is absent, default to 'double'

if nargin == 3
    datatype = 'double';
end

% h5create, h5write won't write empty datasets. Dodge this by only writing
% if array is non-empty. Correspondingly, h5load returns an empty array
% if the dataset doesn't exist.

if ~isempty(array)
        h5create(h5_file, datasetname, size(array), 'datatype', datatype);
        h5write(h5_file, datasetname, array);
end