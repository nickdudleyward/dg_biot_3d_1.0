function rhs_fieldX_e = biot_rhs_3d_e(param, adjoint_calc, dfieldX_dx_e, dfieldX_dy_e, dfieldX_dz_e, face_flux_e)

% Modified for forward adjoint case 2020/07/08 SPE/NFDW
% Modified for pure elastic model 2020/08/06 SPE/NFDW

% Return empty array if there are no elastic elements

if isempty(dfieldX_dx_e)
    rhs_fieldX_e = [];
    return
end

Globals3D;

% Select adjoint or non-adjoint matrices for sparse multiplication loop
% below

if adjoint_calc
    QiA = param.QiA_e_adj;
    QiB = param.QiB_e_adj;
    QiC = param.QiC_e_adj;
else
    QiA = param.QiA_e;
    QiB = param.QiB_e;
    QiC = param.QiC_e;
end

% Compute QiA*dfieldX_dx + QiB*dfieldX_dy + QiC*dfieldX_dz. These are 
% sparse so only add where the entry is non-zero.

rhs_fieldX_e = zeros(Np,numel(param.elastic_elt),param.Nfields);

% Compute QiA*dfieldX_dx + QiB*dfieldX_dy + QiC*dfieldX_dz on each element node.
% These are sparse so only add where the entry is non-zero. This assumes
% that the sparseness pattern is the same on all elements (third index)

% Note on use of squeeze:
% QiA etc are stored as three-dimensional arrays, with element number in
% the third index. dfieldX_dx_e etc. are stored as two-dimensional arrays,
% with element number as column. We want to pointwise multiply each row of
% dfieldX_dx_e by the values in the slice QiA(fld1, fld2, :). For this,
% we need the values in QiA(fld1, fld2, :) to be a row vector (Matlab
% will then automatically multiply all rows), but it has dimensions
% [1, 1, num_elastic_elt]. Squeezing it gives a column vector, so we need
% to squeeze and transpose to give a row vector.

for fld1 = 1:param.Nfields
    volu = zeros(Np,numel(param.elastic_elt));
    for fld2=1:param.Nfields
        if QiA(fld1, fld2, 1) ~= 0
            volu = volu + squeeze(QiA(fld1, fld2, :))'.*dfieldX_dx_e(:,:,fld2);
        end
        if QiB(fld1, fld2, 1) ~= 0
            volu = volu + squeeze(QiB(fld1, fld2, :))'.*dfieldX_dy_e(:,:,fld2);
        end
        if QiC(fld1, fld2, 1) ~= 0
            volu = volu + squeeze(QiC(fld1, fld2, :))'.*dfieldX_dz_e(:,:,fld2);
        end
    end % fld2
    rhs_fieldX_e(:,:,fld1) = -volu - LIFT*(Fscale(:,param.elastic_elt).*face_flux_e(:,:,fld1));
end % fld1
