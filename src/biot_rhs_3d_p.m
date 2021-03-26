function rhs_fieldX_p = biot_rhs_3d_p(param, adjoint_calc, dfieldX_dx_p, dfieldX_dy_p, dfieldX_dz_p, face_flux_p)

% Modified for forward adjoint case 2020/07/08 SPE/NFDW

% Return empty array if there are no poroelastic elements

if isempty(dfieldX_dx_p)
    rhs_fieldX_p = [];
    return
end

Globals3D;

if adjoint_calc
    QiA = param.QiA_p_adj;
    QiB = param.QiB_p_adj;
    QiC = param.QiC_p_adj;
else
    QiA = param.QiA_p;
    QiB = param.QiB_p;
    QiC = param.QiC_p;
end

Globals3D;

rhs_fieldX_p = zeros(Np,numel(param.poroelastic_elt),param.Nfields);

% Compute QiA*dfieldX_dx + QiB*dfieldX_dy + QiC*dfieldX_dz on each element node.
% These are sparse so only add where the entry is non-zero. This assumes
% that the sparseness pattern is the same on all elements (third index)

% Note on use of squeeze:
% QiA_p etc are stored as three-dimensional arrays, with element number in
% the third index. dfieldX_dx_p etc. are stored as two-dimensional arrays,
% with element number as column. We want to pointwise multiply each row of
% dfieldX_dx_p by the values in the slice QiA_p(fld1, fld2, :). For this,
% we need the values in QiA_p(fld1, fld2, :) to be a row vector (Matlab
% will then automatically multiply all rows), but it has dimensions
% [1, 1, num_poroelastic_elt]. Squeezing it gives a column vector, so we #
% need to squeeze and transpose to give a row vector.

%keyboard
for fld1 = 1:param.Nfields
    volu = zeros(Np,numel(param.poroelastic_elt));
    for fld2=1:param.Nfields
        if QiA(fld1, fld2, 1) ~= 0
            volu = volu + squeeze(QiA(fld1, fld2, :))'.*dfieldX_dx_p(:,:,fld2);
        end
        if QiB(fld1, fld2, 1) ~= 0
            volu = volu + squeeze(QiB(fld1, fld2, :))'.*dfieldX_dy_p(:,:,fld2);
        end
        if QiC(fld1, fld2, 1) ~= 0
            volu = volu + squeeze(QiC(fld1, fld2, :))'.*dfieldX_dz_p(:,:,fld2);
        end
    end % fld2
    rhs_fieldX_p(:,:,fld1) = -volu - LIFT*(Fscale(:,param.poroelastic_elt).*face_flux_p(:,:,fld1));    
end % fld1

% dissipative terms are computed in dissipation_lf and dissipation_hf

