function field = combine_memory(field, mem_hf, d, param, zero)

% Pads the memory variables in mem_hf to the full size of the field array 
% then concatenates with the field values in field, to return an augmented
% field array with the memory variables attached at the end.
%
% Both field and mem_hf should be 3-dimensional on input. Concatenated
% structure is reshaped to d dimensions after concatenation, where d is
% 1, 2 or 3.
%
% Intended for saving the augmented field block for communication with the 
% C framework, but could have other uses.
% 
% If the last argument is present and true, use zeros for padding;
% otherwise, use NaN for padding.

% Return field with no change if there are no memory variables.
% If field itself is empty (e.g. called with adjoint field on a 
% non-adjoint run), corresponding mem_hf will surly be empty and this will
% return an empty answer.

if isempty(mem_hf)
    return
end

% Initialise storage with the same number of rows and columns (Np x K) 
% as the main field block, and the number of memory variables (3) in the 
% third dimension.

if nargin == 5 && zero
    mem_inflated = NaN  (size(field,1), size(field, 2), size(mem_hf, 3));
else
    mem_inflated = zeros(size(field,1), size(field, 2), size(mem_hf, 3));
end

% Insert memory variables into initialised storage

mem_inflated(:, param.hf_elt, :) = mem_hf;

% and concatenate the field and memory variables on axis 3

field = cat(3, field, mem_inflated);

% id d is 1 or 2, reshape (without reordering) to d-dimensional array.

switch d
    case 1
        field = field(:);
    case 2
        field = reshape(field, size(field, 1)*size(field, 2), size(field, 3));
    case 3
        % already three dimensional
    otherwise
        error('combine_memory called with dimension %d (must be 1, 2 or 3)', d)
end
     



