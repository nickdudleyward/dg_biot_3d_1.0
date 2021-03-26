function sou = setup_source_moment(param)

Globals3D;

% Find differentiation operators (1 x Np vectors) at the source point.
[Dr_sou,Ds_sou,Dt_sou] = Dmatrices3D(N,param.r_sou, param.s_sou, param.t_sou, V);
% dr/ds etc. on the source element. These are constant
% on columns of the (Np x K) global matrices rx etc.; rx_sou etc are
% scalars
rx_sou = rx(1, param.elm_sou);
ry_sou = ry(1, param.elm_sou);
rz_sou = rz(1, param.elm_sou);
sx_sou = sx(1, param.elm_sou);
sy_sou = sy(1, param.elm_sou);
sz_sou = sz(1, param.elm_sou);
tx_sou = tx(1, param.elm_sou);
ty_sou = ty(1, param.elm_sou);
tz_sou = tz(1, param.elm_sou);
J_sou  =  J(1, param.elm_sou);
% V*V' is inverse mass matrix on reference element; division by 
% J_sou (scalar) gives inverse mass matrix on source element (H&W p.418)
Mi_sou = V*V'/J_sou;

% Assemble source components

Mvec1 = Mi_sou*( (rx_sou*Dr_sou + sx_sou*Ds_sou + tx_sou*Dt_sou)*param.Mxx +...
                 (ry_sou*Dr_sou + sy_sou*Ds_sou + ty_sou*Dt_sou)*param.Mxy +...
                 (rz_sou*Dr_sou + sz_sou*Ds_sou + tz_sou*Dt_sou)*param.Mxz )';
% timo discovered error, Mxz was Mzz before

Mvec2 = Mi_sou*( (rx_sou*Dr_sou + sx_sou*Ds_sou + tx_sou*Dt_sou)*param.Mxy +...
                 (ry_sou*Dr_sou + sy_sou*Ds_sou + ty_sou*Dt_sou)*param.Myy +...
                 (rz_sou*Dr_sou + sz_sou*Ds_sou + tz_sou*Dt_sou)*param.Myz )';
             
             
Mvec3 = Mi_sou*( (rx_sou*Dr_sou + sx_sou*Ds_sou + tx_sou*Dt_sou)*param.Mxz +...
                 (ry_sou*Dr_sou + sy_sou*Ds_sou + ty_sou*Dt_sou)*param.Myz +...
                 (rz_sou*Dr_sou + sz_sou*Ds_sou + tz_sou*Dt_sou)*param.Mzz )';

if param.source_elt_model < param.POROELASTIC
    % the sign is due to definition of the source term -\grad\delta(x-x_sou) in
    % distributional sense
    temp = zeros(param.Nfields, Np);
    temp(param.fld_vx, :) = Mvec1;
    temp(param.fld_vy, :) = Mvec2;
    temp(param.fld_vz, :) = Mvec3;    
    sou  = param.M0 * param.Qi_e(:,:,param.ela_elt_sou) * temp;
elseif param.source_elt_model > param.POROELASTIC
    %keyboard
    temp = -[zeros(7, Np); Mvec1'; Mvec2'; Mvec3'; Mvec1'; Mvec2'; Mvec3'];
    sou  = param.M0 * param.Qi_p(:,:,param.por_elt_sou) * temp;
else
    error('param.source_elt_model = %g not recognised', param.source_elt_model);
end 
             
% sou is (Nfields x Np). This is multiplied by a Ricker wavelet and must 
% be added to all node values of all fields on element elm_sou, i.e. to 
% rhs(:,elm_sou,:) which is (Np x 1 x Nfields), hence transpose and 
% reshape, inserting a degenerate second dimension. 

sou = reshape(sou',Np,1,param.Nfields);

% Note that sou is zero in the first 7 fields since the source acts on the 
% momentum equations and hence only the solid and fluid velocities 
% e.g.
% rhs_field(:,elm_sou,:) = rhs_field(:,elm_sou,:) + wavelet(t)*sou;

