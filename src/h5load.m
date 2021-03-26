function array = h5load(h5file, datasetname)

% h5read insists on absolute names; add a leading / if not
% already present

if datasetname(1) ~= '/'
    datasetname = strcat('/', datasetname);
end

% h5create, h5write won't write empty datasets. h5save dodges this by only 
% writing if array is non-empty; correspondingly, we return an empty array
% if the dataset doesn't exist (or possibly in other error circumstances, 
% too; this catch 'MATLAB:imagesci:h5read:libraryError' could be a pretty 
% blunt instrument!)

try
    array = h5read(h5file, datasetname);
catch ME
    if strcmp(ME.identifier, 'MATLAB:imagesci:h5read:libraryError')
        array = [];
    else
        rethrow(ME);
    end
end
