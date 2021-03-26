function rhs_ker = kernel_rhs(field, rhs_field, field_adj, param)

% 2020/12/17
% Computes right hand side of sensitivity kernels, to be integrated by
% the RK loop in biot_3d
% For the material parameters associated with the mass matrix, the time
% derivatives are computed using the right hand side of the associated field
% ref equations (205)-(207)
% For parameters associated with stiffness matrices, eqns (208)-(209) in
% paper, the spatial derivatives are approximated by their DG
% discretisations

Globals3D;

% initialise storage for output

rhs_ker = zeros(size(field,1), size(field,2), param.num_kernels);

% local names, for convenience

vx_adj = field_adj(:,:,param.fld_vx);
vy_adj = field_adj(:,:,param.fld_vy);
vz_adj = field_adj(:,:,param.fld_vz);
vfx_adj = field_adj(:,:,param.fld_vfx);
vfy_adj = field_adj(:,:,param.fld_vfy);
vfz_adj = field_adj(:,:,param.fld_vfz);

% Density kernels

% this is ker_rho_a. RHS is the dot product of the adjoint solid 
% velocity with the time-derivative of the non-adjoint solid velocity

rhs_ker(:, :, param.fld_ker_rho_a) = ...
    field_adj(:,:,param.fld_vx).*rhs_field(:,:,param.fld_vx) + ...
    field_adj(:,:,param.fld_vy).*rhs_field(:,:,param.fld_vy) + ...
    field_adj(:,:,param.fld_vz).*rhs_field(:,:,param.fld_vz);

% this is ker_rho_f. RHS is the dot product of the adjoint solid 
% velocity with the non-adjoint fluid velocity plus the dot product
% of the adjoint fluid velocity with the non-adjoint solid velocity

rhs_ker(:, :, param.fld_ker_rho_f) = ...
    field_adj(:,:,param.fld_vx) .*rhs_field(:,:,param.fld_vfx) + ...
    field_adj(:,:,param.fld_vy) .*rhs_field(:,:,param.fld_vfy) + ...
    field_adj(:,:,param.fld_vz) .*rhs_field(:,:,param.fld_vfz) + ...
    field_adj(:,:,param.fld_vfx).*rhs_field(:,:,param.fld_vx)  + ...
    field_adj(:,:,param.fld_vfy).*rhs_field(:,:,param.fld_vy)  + ...
    field_adj(:,:,param.fld_vfz).*rhs_field(:,:,param.fld_vz);

% this is ker_m. RHS is the dot product of the adjoint fluid 
% velocity with the time-derivative of the non-adjoint fluid velocity

rhs_ker(:, :, param.fld_ker_m) = ...
    field_adj(:,:,param.fld_vfx).*rhs_field(:,:,param.fld_vfx) + ...
    field_adj(:,:,param.fld_vfy).*rhs_field(:,:,param.fld_vfy) + ...
    field_adj(:,:,param.fld_vfz).*rhs_field(:,:,param.fld_vfz);

% Stiffness kernels

% Differentiate strain tensor in r, s, t  (reference element) coords

d_e11_dr  = Dr*field(:, :, param.fld_e11);
d_e22_dr  = Dr*field(:, :, param.fld_e22);
d_e33_dr  = Dr*field(:, :, param.fld_e33);
d_e12_dr  = Dr*field(:, :, param.fld_e12);
d_e23_dr  = Dr*field(:, :, param.fld_e23);
d_e13_dr  = Dr*field(:, :, param.fld_e13);
d_zeta_dr = Dr*field(:, :, param.fld_zeta);

d_e11_ds  = Ds*field(:, :, param.fld_e11);
d_e22_ds  = Ds*field(:, :, param.fld_e22);
d_e33_ds  = Ds*field(:, :, param.fld_e33);
d_e12_ds  = Ds*field(:, :, param.fld_e12);
d_e23_ds  = Ds*field(:, :, param.fld_e23);
d_e13_ds  = Ds*field(:, :, param.fld_e13);
d_zeta_ds = Ds*field(:, :, param.fld_zeta);

d_e11_dt  = Dt*field(:, :, param.fld_e11);
d_e22_dt  = Dt*field(:, :, param.fld_e22);
d_e33_dt  = Dt*field(:, :, param.fld_e33);
d_e12_dt  = Dt*field(:, :, param.fld_e12);
d_e23_dt  = Dt*field(:, :, param.fld_e23);
d_e13_dt  = Dt*field(:, :, param.fld_e13);
d_zeta_dt = Dt*field(:, :, param.fld_zeta);

d_traceE_dr = d_e11_dr + d_e22_dr + d_e33_dr;
d_traceE_ds = d_e11_ds + d_e22_ds + d_e33_ds;
d_traceE_dt = d_e11_dt + d_e22_dt + d_e33_dt;

rhs_ker(:, :, param.fld_ker_kappa_fr) = - ...
     (vx_adj.*rx + vy_adj.*ry + vz_adj.*rz).*d_traceE_dr - ...
     (vx_adj.*sx + vy_adj.*sy + vz_adj.*sz).*d_traceE_ds - ...
     (vx_adj.*tx + vy_adj.*ty + vz_adj.*tz).*d_traceE_dt;

d_e11_dx = d_e11_dr.*rx + d_e11_ds.*sx + d_e11_dt.*tx;
d_e12_dx = d_e12_dr.*rx + d_e12_ds.*sx + d_e12_dt.*tx;
d_e13_dx = d_e13_dr.*rx + d_e13_ds.*sx + d_e13_dt.*tx;

d_e12_dy = d_e12_dr.*ry + d_e12_ds.*sy + d_e12_dt.*ty;
d_e22_dy = d_e22_dr.*ry + d_e22_ds.*sy + d_e22_dt.*ty;
d_e23_dy = d_e23_dr.*ry + d_e23_ds.*sy + d_e23_dt.*ty;

d_e13_dz = d_e13_dr.*rz + d_e13_ds.*sz + d_e13_dt.*tz;
d_e23_dz = d_e23_dr.*rz + d_e23_ds.*sz + d_e23_dt.*tz;
d_e33_dz = d_e33_dr.*rz + d_e33_ds.*sz + d_e33_dt.*tz;

rhs_ker(:, :, param.fld_ker_mu_fr) = - ...
    vx_adj.*(d_e11_dx + d_e12_dy + d_e13_dz) - ...
    vy_adj.*(d_e12_dx + d_e22_dy + d_e23_dz) - ...
    vz_adj.*(d_e13_dx + d_e23_dy + d_e33_dz) + ...
    rhs_ker(:, :, param.fld_ker_kappa_fr)/3;

% kernels for coupling coefficients

% (210) in paper, intermediate step in computing ker_alpha and ker_M
% below

rhs_ker_alpha_M = ...
     (vx_adj.*rx + vy_adj.*ry + vz_adj.*rz).*d_zeta_dr + ...
     (vx_adj.*sx + vy_adj.*sy + vz_adj.*sz).*d_zeta_ds + ...
     (vx_adj.*tx + vy_adj.*ty + vz_adj.*tz).*d_zeta_dt - ...
     (vfx_adj.*rx + vfy_adj.*ry + vfz_adj.*rz).*d_traceE_dr - ...
     (vfx_adj.*sx + vfy_adj.*sy + vfz_adj.*sz).*d_traceE_ds - ...
     (vfx_adj.*tx + vfy_adj.*ty + vfz_adj.*tz).*d_traceE_dt;

% param.M, param.alpha are K x 1 arrays, holding values per element;
% repmat with transpose replicates these to per-node (Np x K)
 
rhs_ker(:, :, param.fld_ker_alpha) = ...
    repmat(param.M', Np, 1) .* rhs_ker_alpha_M;

rhs_ker(:, :, param.fld_ker_M) = ...
    repmat(param.alpha', Np, 1) .* rhs_ker_alpha_M + ...
    (vfx_adj.*rx + vfy_adj.*ry + vfz_adj.*rz).*d_zeta_dr + ...
    (vfx_adj.*sx + vfy_adj.*sy + vfz_adj.*sz).*d_zeta_ds + ...
    (vfx_adj.*tx + vfy_adj.*ty + vfz_adj.*tz).*d_zeta_dt;

