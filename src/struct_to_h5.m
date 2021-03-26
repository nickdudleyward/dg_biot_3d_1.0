function struct_to_h5(data, h5_file, info, new_name)
% Writes the contents of the struct 'data' to an hdf5 file 'h5_file'
% in the current directory. IF THIS FILE ALREADY EXISTS, IT IS DELETED
% BEFORE WRITING. Field names in the struct map directly to top-level 
% dataset names in the hdf5 file (bur see new_name argument below). 
% Vectors (row or column) are stored as one-dimensional arrays, other 
% arrays are stored in their natural shape. 
%
% Datatype defaults to 'double', but can be overridden by the optional 
% third argument, which should be a struct; if the field in 'data' also 
% exists in 'info', its value in 'info' should be a valid type 
% for h5create, e.g. 'int32', OR the special keyword 'skip' in which case
% the field in 'data' is ignored. 
%
% The dataset name in the HDF5 file defaults to the field name in the
% struct, but can be overridden by the optional fourth argument, which 
% should be a struct; if the field in 'data' also exists in 'new_name',
% the value in 'new-name' is used as the dataset name.
%
% The 'info' and 'new_name' arguments are optional, but it is not possible 
% to specify new_name without info; pass the empty struct struct() as the 
% third argument if required.
%
% Non-numeric data cannot be written to the HDF5 file; such fields are 
% ignored except that, unless they have been marked as 'skip' in 'info', a 
% warning message is printed. Empty matrices cannot (by Matlab) be written 
% to the HDF5 file; such fields are ignored except that, unless they have 
% been marked as 'skip' in 'info', a warning message is printed.

h5_file = absname(h5_file);

if exist(h5_file, 'file') == 2
    delete(h5_file);
    fid = fopen('struct_to_h5.log', 'w');
    fprintf(fid, 'Deleted %s %s\n', pwd, h5_file);
    fclose(fid);
else
    fid = fopen('struct_to_h5.log', 'w');
    fprintf(fid, 'No file %s %s\n', pwd, h5_file);    
    fclose(fid);
end

fields = fieldnames(data);

for f=1:numel(fields)
    
    fname = fields{f};
        
    if nargin >= 3 && isfield(info, fname)
        finfo = info.(fname);
    else
        finfo = '';
    end
    
    if  strcmp(finfo, 'skip')
        continue
    end

    fdata = data.(fname);

    if islogical(fdata)
        fdata = int32(fdata);
        if strcmp(finfo, '')
            finfo = 'int32';
        end
    end
    
    if ~isnumeric(fdata)
        warning('Skipping non-numeric field %s', fname);
        continue
    end
    
    if isempty(fdata)
        warning('Skipping empty field %s', fname);
        continue
    end

    if strcmp(finfo, '')
        ftype = 'double';
    else
        ftype = finfo;
    end
    
    if isvector(fdata) % also true for scalars
        dsize = numel(fdata);
    else
        dsize = size(fdata);
    end
    
    if nargin == 4 && isfield(new_name, fname)
        dname = strcat('/', new_name.(fname));
    else
        dname = strcat('/', fname);
    end

    fprintf('%s %s\n',dname, mat2str(dsize));
    
    if isreal(fdata)
        write_data(h5_file, dname, dsize, ftype, fdata);
    else
        write_data(h5_file, strcat(dname, '_re'), dsize, ftype, real(fdata));
        write_data(h5_file, strcat(dname, '_im'), dsize, ftype, imag(fdata));
    end        
end
end

function write_data(h5_file, dname, dsize, ftype, fdata)
    fid = fopen('struct_to_h5.log', 'a');
    fprintf(fid, '%s %s Field %s\n', pwd, h5_file, dname);
    fclose(fid);
    h5create(h5_file, dname, dsize, 'Datatype', ftype);
    h5write( h5_file, dname, fdata);
end
