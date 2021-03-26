function param = setup_receiver(param)

Globals3D;

% If there are no receiver locations, return empty matrices for the 
% receiver elements and interpolation vectors

if isempty(param.g_rec)
    param.elm_rec = [];
    param.interp_rec = [];
    return
end

% find elements containing receivers and barycentric coordinates
% of receivers within elements
% param.g_rec is an (n x 3) matrix whose rows contain receiver coords

[param.elm_rec, bary_rec] = tsearchn([VX',VY',VZ'], EToV, param.g_rec);

% param.elm_rec is a column vector of n element numbers and bary_rec is an
% (n x 4) matrix of barycentric coordinates. Any coordinate not found in 
% mesh generates a NaN in param.elm_rec

failed_rec = isnan(param.elm_rec);
if any(failed_rec)
    error('Receiver coordinates %s not located in mesh', mat2str(param.g_rec(failed_rec, :)));
end

% Convert barycentrics to r, s, t coordinate system: H&W p.409
% rst_rec is an (n x 3) matrix whose rows contains local element 
% coordinates of receivers

rst_rec = bary_rec*[[-1, -1, -1]; [1, -1, -1]; [-1, 1, -1]; [-1, -1, 1]];

% Vandermonde matrix for evaluating field at each receiver point
% (from field modal values on each receiver element). V_rec is (n x Np)
% for n receivers, i.e. each row holds the interpolation vector for a
% receiver
% Note that Vandermonde3D behaves very strangely if it gets rows instead
% of columns for its r, s, t arguments! Save V_rec in param so 
% setup_adjoint_wavefield_source can use it.
param.V_rec = Vandermonde3D(N, rst_rec(:,1), rst_rec(:,2), rst_rec(:,3));
% Interpolation matrix for evaluating field at each receiver point
% (from field nodal values on each receiver element): global invV maps 
% from nodal to modal values.
param.interp_rec = param.V_rec*invV;

% Find element model type (elastic or poroelastic) for each receiver
% (used in adjoint source)

param.rec_elt_model = param.elt_to_model_type(param.elm_rec);