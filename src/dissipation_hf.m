% Apply high-frequency dissipation. On calling, field_hf and rhs_field_hf 
% should contain all field values on all element nodes of all 
% high-frequency poroelastic elements and the RHS values calculated without
% attenuation; mem_hf should contain the memory variables. On return, 
% rhs_field_hf contains the dissipated values; it should replace the values 
% passed as rhs_field_hf

% Note on use of squeeze:
% QiD_hf is stored as a three-dimensional array, with element number in
% the third index. field_hf etc. are stored as two-dimensional arrays,
% with element number as column. We want to pointwise multiply each row of
% field_hf by the values in the slice QiD_hf(fld1, fld2, :). For this,
% we need the values in QiD_hf(fld1, fld2, :) to be a row vector (Matlab
% will then automatically multiply all rows), but it has dimensions
% [1, 1, num_hf_elt]. Squeezing it gives a column vector, so we need
% to squeeze and transpose to give a row vector.

function [rhs_field_hf, rhs_mem_hf] = dissipation_hf(field_hf, rhs_field_hf, mem_hf, param)

rhs_mem_hf = zeros(size(mem_hf));

% If there are no HF elements, return the unchanged rhs_field_hf and empty
% (as initialised above) ths_mem_hf

if isempty(mem_hf)
    return
end

rhs_field_hf(:,:,8)  = rhs_field_hf(:,:,8)  + squeeze(param.QiD_hf( 8,11,:))'.*field_hf(:,:,11) + squeeze(param.QiD_hf( 8,14,:))'.*mem_hf(:,:,1);
rhs_field_hf(:,:,9)  = rhs_field_hf(:,:,9)  + squeeze(param.QiD_hf( 9,12,:))'.*field_hf(:,:,12) + squeeze(param.QiD_hf( 9,15,:))'.*mem_hf(:,:,2);
rhs_field_hf(:,:,10) = rhs_field_hf(:,:,10) + squeeze(param.QiD_hf(10,13,:))'.*field_hf(:,:,13) + squeeze(param.QiD_hf(10,16,:))'.*mem_hf(:,:,3);

rhs_field_hf(:,:,11) = rhs_field_hf(:,:,11) + squeeze(param.QiD_hf(11,11,:))'.*field_hf(:,:,11) + squeeze(param.QiD_hf(11,14,:))'.*mem_hf(:,:,1);
rhs_field_hf(:,:,12) = rhs_field_hf(:,:,12) + squeeze(param.QiD_hf(12,12,:))'.*field_hf(:,:,12) + squeeze(param.QiD_hf(12,15,:))'.*mem_hf(:,:,2);
rhs_field_hf(:,:,13) = rhs_field_hf(:,:,13) + squeeze(param.QiD_hf(13,13,:))'.*field_hf(:,:,13) + squeeze(param.QiD_hf(13,16,:))'.*mem_hf(:,:,3);

% It is vital to de-couple the memory variables from the physical
% system and solve as ODEs

rhs_mem_hf(:,:,1) =  (param.taue_hf./param.taus_hf - 1)'.*rhs_field_hf(:,:,11) - mem_hf(:,:,1)./param.taus_hf';
rhs_mem_hf(:,:,2) =  (param.taue_hf./param.taus_hf - 1)'.*rhs_field_hf(:,:,12) - mem_hf(:,:,2)./param.taus_hf';
rhs_mem_hf(:,:,3) =  (param.taue_hf./param.taus_hf - 1)'.*rhs_field_hf(:,:,13) - mem_hf(:,:,3)./param.taus_hf';
