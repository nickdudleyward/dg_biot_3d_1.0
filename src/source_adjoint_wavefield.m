function sou_adj = source_adjoint_wavefield(param, rec_num, tstep, s)

% rec_num indexes over receiver number

if isnan(s)
    adj_mag = squeeze(param.source_magnitude_adj_wavefield(rec_num, tstep, :))';
else
    adj_mag_1 = squeeze(param.source_magnitude_adj_wavefield(rec_num, tstep, :))';
    adj_mag_2 = squeeze(param.source_magnitude_adj_wavefield(rec_num, tstep+1, :))';
    adj_mag = adj_mag_1 + s*(adj_mag_2-adj_mag_1);
end

source_components = [zeros(1,7)  adj_mag zeros(1,3)]';

if param.rec_elt_model(rec_num) < param.POROELASTIC
    % TODO check this
    %sou_adj = param.Qi_e(:,:,param.ela_elt_sou)*source_components*param.adjoint_source_IMM(rec_num, :);
elseif param.rec_elt_model(rec_num) > param.POROELASTIC
    % NOTE THAT IN THIS IMPLEMENTATION THE ADJOINT SOURCE IS SIMPLY
    % EVALUATED AT THE RIGHT HAND POINT THROUGHOUT THE RUNGE KUTTA LOOP
    sou_adj = param.Qi_p(:,:,param.elm_rec(rec_num))*source_components*param.adjoint_source_IMM(rec_num, :);
else
    error('param.source_elt_model = %g not recognised', param.source_elt_model);
end    

% sou_adj is (Nfields x Np). This must be added to all node values of all 
% fields on element elm_rec(rec_num), i.e. to 
% rhs(:,elm_rec(rec_num),:) which is (Np x 1 x Nfields), hence transpose and 
% reshape, inserting a degenerate second dimension. 

% (avoiding reference to global Np as loading Globals3D causes warning
% messages)

sou_adj = reshape(sou_adj', [], 1, param.Nfields);