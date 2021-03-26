function sou = setup_source_force(param)
% This models surface point force param.force = [F_1 F_2 F_3]

Globals3D;

% V_sou is P(x_sou), the Lagrange poloynomials evaluated at source
% location, see H&W, page 47 ff, source_components is (Nfields x1)
V_sou = Vandermonde3D(N, param.r_sou, param.s_sou, param.t_sou)';
tmp = 1/J(1, param.elm_sou)*V*V_sou;

if param.source_elt_model < param.POROELASTIC
    source_components = zeros(param.Nfields, 1);
    source_components(param.fld_vx:param.fld_vz) = param.force;
    sou = param.Qi_e(:,:,param.ela_elt_sou)*source_components*tmp';
elseif param.source_elt_model > param.POROELASTIC
    source_components = [zeros(1,7) param.force param.force]';
    sou = param.Qi_p(:,:,param.por_elt_sou)*source_components*tmp';  
else
    error('param.source_elt_model = %g not recognised', param.source_elt_model);
end    

% sou is (Nfields x Np). This is multiplied by a Ricker wavelet and must 
% be added to all node values of all fields on element elm_sou, i.e. to 
% rhs(:,elm_sou,:) which is (Np x 1 x Nfields), hence transpose and 
% reshape, inserting a degenerate second dimension. 

sou = reshape(sou',Np,1,param.Nfields);


 