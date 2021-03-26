function [rhs_fieldX, rhs_mem_hfX] = biot_rhs_3d(RKtime, adjoint_calc, fieldX, mem_hfX, param)

% Modified for forward adjoint case 2020/07/08 SPE/NFDW
% fieldX etc can be wavefield (adjoint_calc false) or adjoint wavefield
% (adjoint_calc true)

persistent tstep;
if isempty(tstep)
    tstep = 1;
end

Globals3D;

persistent call_number;
if isempty(call_number)
    call_number = 1;
else
    call_number = call_number + 1;
end

% Start off with the interior and exterior values, without boundary values
% so no field differences at boundary nodes

% size(fieldX) at this point is [Np, K, param.Nfields]; reshape so we can
% use global element node numbers to extract interior and exterior values

fieldX = reshape(fieldX, K*Np, param.Nfields);

fieldXM = fieldX(vmapM, :);
fieldXP = fieldX(vmapP, :);

% change field back to original shape, as needed later

fieldX = reshape(fieldX, Np, K, param.Nfields);

% apply boundary value multipliers to fieldXP
% F denotes far/outflow boundary condition
% W denotes wall or free surface bc
% need intersection when domains are separated

bmultF = -ones(1,param.Nfields);
bmultW = [-ones(1,param.num_str_fields), ones(1,param.num_vel_fields)];

% this fudge makes upper face absorbing for checking
% bmultW = [-ones(1,param.num_str_fields), -ones(1,param.num_vel_fields)];

fieldXP(mapF,:) = fieldXP(mapF,:).*repmat(bmultF, numel(mapF), 1);
fieldXP(mapW,:) = fieldXP(mapW,:).*repmat(bmultW, numel(mapW), 1);

% Apply Dirichlet BCs, used in convergence tests
% Need isempty test because param.DirichletBC is probably undefined if
% there are no Dirichlet BCs

if ~isempty(mapD)
    fieldXP(mapD, :) = param.DirichletBC(RKtime, x(vmapD), y(vmapD), z(vmapD), param);
end

% Now find fluxes on faces

% get_face_flux_* returns empty if there are no faces of the relevant type,
% so no change is made to face_flux in that case

% Get face fluxes for all three interface types (elastic-elastic, 
% poroelastic-poroelastic and elastic-poroelastic) and insert into
% face_flux to build a global face flux matrix

face_flux = NaN(Nfaces*Nfp*K, param.Nfields);

if adjoint_calc
    face_flux_ee = get_face_flux_ee_adj(fieldXM, fieldXP, param);
    face_flux_pp = get_face_flux_pp_adj(fieldXM, fieldXP, param);
else
    face_flux_ee = get_face_flux_ee(fieldXM, fieldXP, param);
    face_flux_pp = get_face_flux_pp(fieldXM, fieldXP, param);
end

face_flux(param.ee_face_node,:) = face_flux_ee;
face_flux(param.pp_face_node,:) = face_flux_pp;

% There are potentially two ways to calculate face_flux_ep. If the nodes 
% are ordered in the right way, can use get_face_flux_ep_pe or 
% get_face_flux_ep. In this case, do both and compare answers as a test.
% If the nodes are not ordered in this way, get_face_flux_ep_pe will only
% return the pe face flux and get_face_flux_ep will be needed to get the
% ep flux.

if param.align_ep_pe_nodes
    if adjoint_calc
        [face_flux_ep, face_flux_pe] = get_face_flux_ep_pe_adj(fieldXM, fieldXP, param);
        face_flux_ep_new = get_face_flux_ep_adj(fieldXM, fieldXP, param);
    else
        [face_flux_ep, face_flux_pe] = get_face_flux_ep_pe(fieldXM, fieldXP, param);
        face_flux_ep_new = get_face_flux_ep(fieldXM, fieldXP, param);
    end
    face_flux(param.pe_face_node,:) = face_flux_pe;
    face_flux(param.ep_face_node,:) = face_flux_ep;
    fprintf('flux difference %.5e / %.5e\n', ...
        max(abs(face_flux_ep(:)-face_flux_ep_new(:))), ...
        max(abs(face_flux_ep(:))));
else
    if adjoint_calc    
        [~, face_flux_pe] = get_face_flux_ep_pe_adj(fieldXM, fieldXP, param);
        face_flux_ep = get_face_flux_ep_adj(fieldXM, fieldXP, param);
    else
        [~, face_flux_pe] = get_face_flux_ep_pe(fieldXM, fieldXP, param);
        face_flux_ep = get_face_flux_ep(fieldXM, fieldXP, param);
    end
    face_flux(param.pe_face_node,:) = face_flux_pe;
    face_flux(param.ep_face_node,:) = face_flux_ep;
end

% This is just a quick check that every entry in face_flux has been
% assigned to

if any(isnan(face_flux(:)))
    error('Some face_flux entries unassigned')
end

% Convert from 2D array indexed by global face node number and field to 
% 3D array indexed by local face node number, element number and 
% field number.

face_flux = reshape(face_flux, Nfaces*Nfp, K, param.Nfields);

% Now look at fields on the element body

% Differentiate all fields

dfieldX_dx = zeros(size(fieldX));
dfieldX_dy = zeros(size(fieldX));
dfieldX_dz = zeros(size(fieldX));

% TODO: stitch this together as one vector operation, instead of
% looping over fields

for fld=1:param.Nfields
    % compute local derivatives
    dfieldX_dr = Dr*fieldX(:,:,fld);
    dfieldX_ds = Ds*fieldX(:,:,fld);
    dfieldX_dt = Dt*fieldX(:,:,fld);
    % piece local derivatives together to compute physical
    % derivatives
    dfieldX_dx(:,:,fld) = rx.*dfieldX_dr + sx.*dfieldX_ds + tx.*dfieldX_dt;
    dfieldX_dy(:,:,fld) = ry.*dfieldX_dr + sy.*dfieldX_ds + ty.*dfieldX_dt;
    dfieldX_dz(:,:,fld) = rz.*dfieldX_dr + sz.*dfieldX_ds + tz.*dfieldX_dt;
end

% Compute RHS (excluding source) on both element types (elastic and
% poroelastic) and insert into rhs_fieldX to build global RHS
% (excluding source)

rhs_fieldX = NaN(size(fieldX));

% biot_rhs_3d_* returns empty if there are no elements of the relevant type,
% so no change is made to rhs_fieldX in that case

rhs_fieldX(:,param.elastic_elt,:) = biot_rhs_3d_e(param, adjoint_calc, ...
    dfieldX_dx(:,param.elastic_elt,:), dfieldX_dy(:,param.elastic_elt,:), ...
    dfieldX_dz(:,param.elastic_elt,:), face_flux(:,param.elastic_elt,:));

rhs_fieldX(:,param.poroelastic_elt,:) = biot_rhs_3d_p(param, adjoint_calc, ...
    dfieldX_dx(:,param.poroelastic_elt,:), dfieldX_dy(:,param.poroelastic_elt,:), ...
    dfieldX_dz(:,param.poroelastic_elt,:), face_flux(:,param.poroelastic_elt,:));

% This is just a quick check that every entry in rhs_field has been
% assigned to

if any(isnan(rhs_fieldX(:)))
    error('Some rhs_fieldX entries unassigned')
end

% If we have any low-frequency elements, apply low-frequency dissipation to
% them, unless we're using operator splitting or IMEX to handle this, in
% which case param.include_lf will be set false

if param.include_lf && ~isempty(param.lf_elt)
    %fprintf('RHS LF dissipation\n');
    rhs_fieldX_lf = dissipation_lf(fieldX(:,param.lf_elt,:), rhs_fieldX(:,param.lf_elt,:), param);
    rhs_fieldX(:,param.lf_elt,:) = rhs_fieldX_lf;
end

% Apply high-frequency dissipation to any HF elements and find new memory 
% variable values. Memory variables remain empty if theer are no HF
% elements.

[rhs_fieldX_hf, rhs_mem_hfX] = dissipation_hf(fieldX(:,param.hf_elt,:), rhs_fieldX(:,param.hf_elt,:), mem_hfX, param);
rhs_fieldX(:,param.hf_elt,:) = rhs_fieldX_hf;

% Add source, if present and required (i.e. not if we're computing the 
% adjoint wavefield for the adjoint method)

if ~(param.run_adjoint_method && adjoint_calc) && ~isempty(param.elm_sou)
    rhs_fieldX(:,param.elm_sou,:) = rhs_fieldX(:,param.elm_sou,:) +...
        ricker_wavelet(RKtime, param.t0_sou, param.freq_sou)*param.sou;
end

% source for adjoint wavefield, if we're computing the adjoint wavefield
% for the adjoint method)

if param.run_adjoint_method && adjoint_calc && ~isempty(param.elm_rec)
    
    % abbreviate for convenience
    step_time = param.forward_step_time;

    % Find the step times from the previous run bracketing RKtime.
    % First check that RKtime is within the previous run time.
    
    if RKtime < step_time(1) || RKtime > step_time(end)
        error('RKtime %g out of forward step time range %g ... %g\n', RKtime, ...
            step_time(1), step_time(end))
    end

    % Detect if we're at the right-hand endpoint of step_time
    
    if RKtime == step_time(end)
        tstep = numel(step_time);
       
    % Otherwise, see if the previous value of tstep is still correct
    % (stages within one RK step will have the same tstep)
    % (tstep is initialised to 1 so step_time(tstep+1) is well-defined on
    % first run)
    
    elseif tstep == numel(step_time) || RKtime < step_time(tstep) || RKtime >= step_time(tstep+1)

    % Otherwise, find the last index of a step time <= RKtime, so
    % step_time(tstep) <= RKtime < step_time(tstep+1) (unless tstep is last index)

        tstep = find(step_time <= RKtime, 1, 'last');
    end

    % location of RKtime in the interval step_time(tstep)..step_time(tstep+1)
    % or NaN if we're at the righthand endpoint
    
    if tstep == numel(step_time)
        s = NaN;
    else
        s = (RKtime - step_time(tstep))/(step_time(tstep+1)-step_time(tstep));
    end
    for ii=1:param.num_rec
        source_adj = source_adjoint_wavefield(param, ii, tstep, s);
        rhs_fieldX(:,param.elm_rec(ii),:) = rhs_fieldX(:,param.elm_rec(ii),:) +...
            source_adj;       
    end
end

%fprintf('Reached end of biot_rhs_3d, step number %d\n', call_number);
%keyboard
