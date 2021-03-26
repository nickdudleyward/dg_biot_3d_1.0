% Apply low-frequency dissipation. On calling, field_lf and rhs_field_lf 
% should contain all field values on all element nodes of all low-frequency 
% poroelastic elements and RHS values compted without attenuation. On 
% return rhs_field_lf contains the dissipated values; it should replace 
% the values passed as rhs_field_lf

function rhs_field_lf = dissipation_lf(field_lf, rhs_field_lf, param)

%fprintf('Explicit LF dissipation\n');

% Note on use of squeeze:
% Qi_p is stored as a three-dimensional array, with element number in
% the third index. field_lf(:,:,fld) is a two-dimensional array,
% with element number as column. We want to pointwise multiply each row of
% field_lf(:,:,fld) by the values in the slice Qi_p(fld1, fld2, :). For this,
% we need the values in Qi_p(fld1, fld2, :) to be a row vector (Matlab
% will then automatically multiply all rows), but it has dimensions
% [1, 1, num_lf_elt]. Squeezing it gives a column vector, so we
% need to squeeze and transpose to give a row vector.
%keyboard
% update solid velocity parameters
for fld = 8:10
    rhs_field_lf(:,:,fld) = rhs_field_lf(:,:,fld) -...
        (squeeze(param.Qi_lf(fld, fld+3, :)).*param.eta_lf./param.k_lf)'.*field_lf(:,:,fld+3);
end
% update fluid velocity parameters
for fld = 11:13
    rhs_field_lf(:,:,fld) = rhs_field_lf(:,:,fld) -...
        (squeeze(param.Qi_lf(fld, fld  , :)).*param.eta_lf./param.k_lf)'.*field_lf(:,:,fld);
end
