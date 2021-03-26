function param = build_param(param)
% Builds data structure with material parameters for each element for
% elastic and poroelastic cases

Globals3D;

% overwritten by save and restart code

param.offset_time = 0;

% For convergence tests, plane_wave_*.m create plane waves propagating
% in this direction

param.pwx = 1;
param.pwy = 1;
param.pwz = 1;

% Make this not parallel or normal to any faces in the triv_mesh output, 
% in case there's any kind of interaction there

param.pwx = 0.9;
param.pwy = 1.0;
param.pwz = 1.1;

% How many steps between printing diagnostic messages. Do not overwrite
% if already set.
% TODO this is now set in driver_test_new; once this is established,
% can be removed from here.

if ~isfield(param, 'steps_per_report')
    param.steps_per_report = 10;
end

% These are just for convenience

num_elastic_elt     = numel(param.elastic_elt);
num_poroelastic_elt = numel(param.poroelastic_elt);
num_hf_elt          = numel(param.hf_elt);
num_lf_elt          = numel(param.lf_elt);
num_inviscid_elt    = numel(param.inviscid_elt);

num_ee_face_node    = numel(param.ee_face_node);
num_pp_face_node    = numel(param.pp_face_node);
num_pe_face_node    = numel(param.pe_face_node);

% THIS IS SPECIFIC TO THE DEFINITION OF PHYSICAL PROPERTIES BY DOMAIN AND
% MIGHT BE BETTER AS A SEPARATE FUNCTION IN FUTURE

% For each field in struct array phys, allocate storage: a row vector, 
% indexed by element number. Initialise this with NaNs, in the hope that
% this will make errors in data allocation obvious.

% Note that field names apply to the whole struct array, not to individual
% entries; if param.phys(dom).name has not been assigned, it holds an empty
% array.

names = fieldnames(param.phys);
for fn = 1:numel(names)
    param.(names{fn}) = NaN(K, 1);
end

% Each entry in struct array param.phys corresponds to a domain. For each field,
% replicate the scalar value in param.phys across all elements of the
% corresponding domain

for dom = 1:numel(param.phys)
    dom_mask = (param.elt_to_domain == dom);
    for fn = 1:numel(names)
        if ~isempty(param.phys(dom).(names{fn}))
            param.(names{fn})(dom_mask) = param.phys(dom).(names{fn});
        end
    end
end

% END SPECIFIC TO THE DEFINITION OF PHYSICAL PROPERTIES BY DOMAIN

% derived poroelastic parameters. These need to be done here before the
% element and face node values are extracted. 

if num_poroelastic_elt > 0
    param.alpha  = 1-param.k_fr./param.k_s;
    param.M      = param.k_s./((1-param.k_fr./param.k_s)-param.phi.*(1-param.k_s./param.k_f));
    param.lambda = param.k_fr+param.alpha.^2.*param.M-2/3*param.mu_fr;
    param.m      = param.rho_f.*param.tau./param.phi;
    param.rho_a  = (1-param.phi).*param.rho_s + param.phi.*param.rho_f;
end

% from the above, we can compute wavespeeds needed for the 'd' 
% coefficients to build the beta's for the P-wave flux terms
% see equations 106-109. Again, this is done here before the poroelastic 
% element and face node values are extracted

if num_poroelastic_elt > 0
    param = get_inviscid_wavespeeds(param);
end

% SWITCH WAVESPEEDS TO GET DOWNWIND SCHEME FOR TIME REVERSAL
% (poroelastic only; elastic not yet implemented)

if param.num_steps < 0 || (param.num_steps == 0 && param.time_width < 0)
    param.cpI_i  = -param.cpI_i;
    param.cpII_i = -param.cpII_i;
    param.cs_i   = -param.cs_i;
end

% All properties are now defined on all elements (as NaN if they
% don't apply to that element). For elastic and poroelastic properties,
% extract values only on relevant elements and (for flux calculations) 
% face nodes.

% Do face nodes first. vmapP and vmapM contain element node numbers;
% emapM and emapP contain corresponding element numbers. This allows us to
% easily inflate data indexed by element number to data indexed by face 
% node number. 

emapM = floor((vmapM - 1) / Np) + 1;
emapP = floor((vmapP - 1) / Np) + 1;

% Elastic properties on elastic-elastic face nodes

if num_ee_face_node > 0

    emap_eeM = emapM(param.ee_face_node);
    emap_eeP = emapP(param.ee_face_node);

    param.rho_eeM    = param.rho_ela   (emap_eeM);
    param.cp_eeM     = param.cp_ela    (emap_eeM);
    param.cs_eeM     = param.cs_ela    (emap_eeM);
    param.mu_eeM     = param.mu_ela    (emap_eeM);
    param.lambda_eeM = param.lambda_ela(emap_eeM);

    param.rho_eeP    = param.rho_ela   (emap_eeP);
    param.cp_eeP     = param.cp_ela    (emap_eeP);
    param.cs_eeP     = param.cs_ela    (emap_eeP);
    param.mu_eeP     = param.mu_ela    (emap_eeP);
    param.lambda_eeP = param.lambda_ela(emap_eeP);

end % if num_ee_face_node > 0

% Poroelastic properties on poroelastic-poroelastic face nodes

if num_pp_face_node > 0

    emap_ppM = emapM(param.pp_face_node);
    emap_ppP = emapP(param.pp_face_node);

    param.alpha_ppM  = param.alpha(emap_ppM);
    param.mu_fr_ppM  = param.mu_fr(emap_ppM);
    param.lambda_ppM = param.lambda(emap_ppM);
    param.M_ppM      = param.M(emap_ppM);
    param.m_ppM      = param.m(emap_ppM);
    param.rho_f_ppM  = param.rho_f(emap_ppM);
    param.cs_i_ppM   = param.cs_i(emap_ppM);
    param.cpI_i_ppM  = param.cpI_i(emap_ppM);
    param.cpII_i_ppM = param.cpII_i(emap_ppM);
    param.gamma1_ppM = param.gamma1(emap_ppM);
    param.gamma2_ppM = param.gamma2(emap_ppM);
    param.cpI_i_ppM  = param.cpI_i(emap_ppM);
    param.cpII_i_ppM = param.cpII_i(emap_ppM);

    param.alpha_ppP  = param.alpha(emap_ppP);
    param.mu_fr_ppP  = param.mu_fr(emap_ppP);
    param.lambda_ppP = param.lambda(emap_ppP);
    param.M_ppP      = param.M(emap_ppP);
    param.m_ppP      = param.m(emap_ppP);
    param.rho_f_ppP  = param.rho_f(emap_ppP);
    param.cs_i_ppP   = param.cs_i(emap_ppP);
    param.cpI_i_ppP  = param.cpI_i(emap_ppP);
    param.cpII_i_ppP = param.cpII_i(emap_ppP);
    param.gamma1_ppP = param.gamma1(emap_ppP);
    param.gamma2_ppP = param.gamma2(emap_ppP);
    param.cpI_i_ppP  = param.cpI_i(emap_ppP);
    param.cpII_i_ppP = param.cpII_i(emap_ppP);

end % if num_pp_face_node > 0

% Elastic and poroelastic properties on elastic-poroelastic face nodes
    
if num_pe_face_node > 0

    % Elastic properties on elastic interior and poroelastic exterior face 
    % nodes
    % These are all interior, i.e. M
    
    emap_peM = emapM(param.pe_face_node);

    param.mu_peM     = param.mu_ela(emap_peM);
    param.lambda_peM = param.lambda_ela(emap_peM);
    param.cp_peM     = param.cp_ela(emap_peM);
    param.cs_peM     = param.cs_ela(emap_peM);
    
    % Poroelastic properties on elastic interior and poroelastic exterior 
    % face nodes
    % These are all exterior, i.e. P

    emap_peP = emapP(param.pe_face_node);

    param.mu_fr_peP  = param.mu_fr(emap_peP);
    param.lambda_peP = param.lambda(emap_peP);
    param.alpha_peP  = param.alpha(emap_peP);
    param.rho_f_peP  = param.rho_f(emap_peP);
    param.M_peP      = param.M(emap_peP);
    param.m_peP      = param.m(emap_peP);
    param.gamma1_peP = param.gamma1(emap_peP);
    param.gamma2_peP = param.gamma2(emap_peP);
    param.cs_i_peP   = param.cs_i(emap_peP);
    param.cpI_i_peP  = param.cpI_i(emap_peP);
    param.cpII_i_peP = param.cpII_i(emap_peP);
    
    % BEGIN NICK
    % THIS IS FOR get_face_flux_ep WITH PORO ELEMENTS ON INTERIOR
    % Poroelastic properties on poroelastic interior and elastic exterior 
    % face nodes
    % These are now interior, ie. M
    
    emap_epM = emapM(param.ep_face_node);
    
    param.mu_fr_epM  = param.mu_fr(emap_epM);
    param.lambda_epM = param.lambda(emap_epM);
    param.alpha_epM  = param.alpha(emap_epM);
    param.rho_f_epM  = param.rho_f(emap_epM);
    param.M_epM      = param.M(emap_epM);
    param.m_epM      = param.m(emap_epM);
    param.gamma1_epM = param.gamma1(emap_epM);
    param.gamma2_epM = param.gamma2(emap_epM);
    param.cs_i_epM   = param.cs_i(emap_epM);
    param.cpI_i_epM  = param.cpI_i(emap_epM);
    param.cpII_i_epM = param.cpII_i(emap_epM);
    
    % Elastic properties on poroelastic interrior and elastic exterior face 
    % nodes
    % These are all exterior, i.e. P
    
    emap_epP = emapP(param.ep_face_node);
    
    param.mu_epP     = param.mu_ela(emap_epP);
    param.lambda_epP = param.lambda_ela(emap_epP);
    param.cp_epP     = param.cp_ela(emap_epP);
    param.cs_epP     = param.cs_ela(emap_epP);
    % END NICK

end % end % if num_pe_face_node > 0


% Elastic properties on elastic elements (do this after face properties
% because it overwrites original data)

if num_elastic_elt > 0

    param.rho_ela    = param.rho_ela   (param.elastic_elt);
    param.cp_ela     = param.cp_ela    (param.elastic_elt);
    param.cs_ela     = param.cs_ela    (param.elastic_elt);
    param.mu_ela     = param.mu_ela    (param.elastic_elt);
    param.lambda_ela = param.lambda_ela(param.elastic_elt);

end % if num_elastic_elt > 0

% Poroelastic properties on poroelastic elements (do this after face 
% properties because it overwrites original data)

if num_poroelastic_elt > 0

    param.alpha  = param.alpha (param.poroelastic_elt);
    param.mu_fr  = param.mu_fr (param.poroelastic_elt);
    param.lambda = param.lambda(param.poroelastic_elt);
    param.M      = param.M     (param.poroelastic_elt);
    param.m      = param.m     (param.poroelastic_elt);
    param.rho_a  = param.rho_a (param.poroelastic_elt);
    param.rho_f  = param.rho_f (param.poroelastic_elt);
    param.eta    = param.eta   (param.poroelastic_elt);

    % Bodged these into conditionals, should move them into
    % appropriate LF/HF only setup below
    
    % These are for LF elements only, for use in dissipation_lf
    if num_lf_elt > 0
        param.eta_lf = param.eta (param.lf_elt);
        param.k_lf   = param.k_lf(param.lf_elt);
    end

    % These are for HF elements only, for use in dissipation_hf
    if num_hf_elt > 0
        param.eta_hf  = param.eta (param.hf_elt);
        param.k_hf    = param.k_hf(param.hf_elt);
        param.taus_hf = param.taus(param.hf_elt);
        param.taue_hf = param.taue(param.hf_elt);
    end

end % if num_poroelastic_elt > 0

% assemble poroelastic elastic stiffness and mass matrices on each 
% poroelastic element. These are Nfields square matrices; the third index 
% in the array is the element number within the list of poroelastic 
% elements.

if num_poroelastic_elt > 0

    A_p = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);

    A_p( 8, 1, :) = 2*param.mu_fr+param.lambda;
    A_p(11, 1, :) = param.M.*param.alpha;
    A_p( 8, 2, :) = param.lambda;
    A_p(11, 2, :) = param.M.*param.alpha;
    A_p( 8, 3, :) = param.lambda;
    A_p(11, 3, :) = param.M.*param.alpha;
    A_p( 9, 4, :) = 2*param.mu_fr;
    A_p(10, 6, :) = 2*param.mu_fr;
    A_p( 8, 7, :) = -param.alpha.*param.M;
    A_p(11, 7, :) = -param.M;
    A_p( 1, 8, :) = 1;
    A_p( 4, 9, :) = 0.5;
    A_p( 6,10, :) = 0.5;
    A_p( 7,11, :) = -1;
    A_p = -A_p; % this is to make consistent with the hyp system formulation d/dt + d/dx(aq) +... =  sources
    param.A_p = A_p;

    B_p = zeros(param.Nfields,param.Nfields, num_poroelastic_elt);
    B_p( 9, 1, :) = param.lambda;
    B_p(12, 1, :) = param.M.*param.alpha;
    B_p( 9, 2, :) = 2*param.mu_fr+param.lambda;
    B_p(12, 2, :) = param.M.*param.alpha;
    B_p( 9, 3, :) = param.lambda;
    B_p(12, 3, :) = param.M.*param.alpha;
    B_p( 8, 4, :) = 2*param.mu_fr;
    B_p(10, 5, :) = 2*param.mu_fr;
    B_p( 9, 7, :) = -param.alpha.*param.M;
    B_p(12, 7, :) = -param.M;
    B_p( 4, 8, :) = 0.5;
    B_p( 2, 9, :) = 1;
    B_p( 5,10, :) = 0.5;
    B_p( 7,12, :) = -1;
    B_p = -B_p;
    param.B_p = B_p;

    C_p = zeros(param.Nfields,param.Nfields, num_poroelastic_elt);
    C_p(10, 1, :) = param.lambda;
    C_p(13, 1, :) = param.M.*param.alpha;    
    C_p(10, 2, :) = param.lambda;
    C_p(13, 2, :) = param.M.*param.alpha;    
    C_p(10, 3, :) = 2*param.mu_fr+param.lambda;
    C_p(13, 3, :) = param.M.*param.alpha;    
    C_p( 9, 5, :) = 2*param.mu_fr;    
    C_p( 8, 6, :) = 2*param.mu_fr;    
    C_p(10, 7, :) = -param.alpha.*param.M;
    C_p(13, 7, :) = -param.M;
    C_p(3, 10, :) = 1;    
    C_p( 5, 9, :) = 0.5;    
    C_p( 6, 8, :) = 0.5;
    C_p( 7,13, :) = -1;
    C_p = -C_p;
    param.C_p = C_p;

    % Diffusivity for low frequency case

    D_p = zeros(param.Nfields, param.Nfields, num_lf_elt);
    if num_lf_elt > 0
        D_p(11,11, :) = -param.eta_lf./param.k_lf;
        D_p(12,12, :) = -param.eta_lf./param.k_lf;
        D_p(13,13, :) = -param.eta_lf./param.k_lf;
    end
    param.D_p = D_p;
    
    % poroelastic mass matrix 
    Q_p = zeros(param.Nfields,param.Nfields, num_poroelastic_elt);
    Q_p (1, 1, :) = 1;
    Q_p( 2, 2, :) = 1;
    Q_p( 3, 3, :) = 1;
    Q_p( 4, 4, :) = 1;
    Q_p( 5, 5, :) = 1;
    Q_p( 6, 6, :) = 1;
    Q_p( 7, 7, :) = 1;
    Q_p( 8, 8, :) = param.rho_a;
    Q_p(11, 8, :) = param.rho_f;    
    Q_p( 9, 9, :) = param.rho_a;
    Q_p(12, 9, :) = param.rho_f;    
    Q_p(10,10, :) = param.rho_a;
    Q_p(13,10, :) = param.rho_f;
    Q_p( 8,11, :) = param.rho_f;
    Q_p(11,11, :) = param.m;
    Q_p( 9,12, :) = param.rho_f;
    Q_p(12,12, :) = param.m;
    Q_p(10,13, :) = param.rho_f;
    Q_p(13,13, :) = param.m;
    param.Q_p = Q_p;

    % Now we need Q^{-1}A, Q^{-1}B, Q^{-1}C for each poroelastic element. 
    % Handle each element individually.
    
    % Also need Q^{-1} itself, for the LF dissipative case and for any
    % poroelastic (not necessarily LF) element containing a source.
    % This is a pain, and is done very inefficiently here. We compute Qi_p
    % for all poroelastic elements and for each LF element take a copy.
    % This still makes finding the Qi_p for the source element a pain
    % because param.elm_sou is an index into all elements but param.Qi_p
    % contains only poroelastic elements

    param.Qi_p   = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);    
    param.QiA_p  = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);
    param.QiB_p  = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);
    param.QiC_p  = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);

    param.Qi_lf  = zeros(param.Nfields, param.Nfields, num_lf_elt);

    elt_lf = 1;
    for elt = 1:num_poroelastic_elt
        param.Qi_p(:,:,elt) = inv(Q_p(:,:,elt));
        param.QiA_p(:,:,elt) =  Q_p(:,:,elt) \ A_p(:,:,elt);
        param.QiB_p(:,:,elt) =  Q_p(:,:,elt) \ B_p(:,:,elt);
        param.QiC_p(:,:,elt) =  Q_p(:,:,elt) \ C_p(:,:,elt);
        % if we have a LF eleemnt, copy Qi to Qi_lf
        if param.elt_to_model_type(elt) == param.LOW_FREQUENCY
            param.Qi_lf(:,:,elt_lf) = param.Qi_p(:,:,elt);
            elt_lf  = elt_lf + 1;
        end
    end % for elt = 1:num_poroelastic_elt
    
end % if num_poroelastic_elt > 0

if num_hf_elt > 0

    % Storage for high-frequency versions, with three extra rows and
    % columns to account for memory variables.
    
    % These are used for analytic plane waves and QiD is used in
    % dissipation_hf
    
    param.A_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
    param.B_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
    param.C_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
    param.D_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
    param.Q_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
    param.QiD_hf = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
    
    % Extend A_p, B_p, C_p to high-frequency versions adding three zero 
    % rows and columns and the bottom and right   
    
    param.A_hf(1:param.Nfields, 1:param.Nfields, :) = A_p(:,:,param.hf_within_poro);
    param.B_hf(1:param.Nfields, 1:param.Nfields, :) = B_p(:,:,param.hf_within_poro);
    param.C_hf(1:param.Nfields, 1:param.Nfields, :) = C_p(:,:,param.hf_within_poro);

    % high frequency mass matrix, extended to include memory variables
    param.Q_hf(1:param.Nfields, 1:param.Nfields, :) = Q_p(:,:,param.hf_within_poro);
    param.Q_hf(14,14,:) = 1;
    param.Q_hf(15,15,:) = 1;
    param.Q_hf(16,16,:) = 1;
    param.Q_hf(14,11,:) = -(param.taue_hf./param.taus_hf-1);
    param.Q_hf(15,12,:) = -(param.taue_hf./param.taus_hf-1);
    param.Q_hf(16,13,:) = -(param.taue_hf./param.taus_hf-1);

    % high frequency diffusivity matrix , extended to include memory variables
    param.D_hf = zeros(param.Nfields, param.Nfields, num_hf_elt);
    param.D_hf(11,11,:)   = -param.eta_hf./param.k_hf;
    param.D_hf(12,12,:)   = -param.eta_hf./param.k_hf;
    param.D_hf(13,13,:)   = -param.eta_hf./param.k_hf;
    param.D_hf(11,14,:)   = -param.eta_hf./param.k_hf;
    param.D_hf(12,15,:)   = -param.eta_hf./param.k_hf;
    param.D_hf(13,16,:)   = -param.eta_hf./param.k_hf;
    param.D_hf(14,14,:)   = -1/param.taus_hf;
    param.D_hf(15,15,:)   = -1/param.taus_hf;
    param.D_hf(16,16,:)   = -1/param.taus_hf;

    for elt=1:num_hf_elt
        param.QiD_hf(:,:,elt) = param.Q_hf(:,:,elt) \ param.D_hf(:,:,elt);
    end
    
end % if num_hf_elt > 0

if num_pp_face_node > 0
    
    %flux stuff, refer eqs (106)-(109) in paper. 

    % this builds unnumbered matrix in \S 3.4. It's not given a name there, but
    % its inverse is d so call it di
    % first column
    
    di_p = zeros(4,4,num_pp_face_node);
    
    di_p(1,1,:) = 2*param.mu_fr_ppM+param.lambda_ppM+...
        param.alpha_ppM.*param.M_ppM.*param.gamma1_ppM;
    di_p(2,1,:) = param.M_ppM.*(param.alpha_ppM+param.gamma1_ppM);
    di_p(3,1,:) = param.cpI_i_ppM;
    di_p(4,1,:) = param.gamma1_ppM.*param.cpI_i_ppM;
    % second column
    di_p(1,2,:) = 2*param.mu_fr_ppM+param.lambda_ppM+...
        param.alpha_ppM.*param.M_ppM.*param.gamma2_ppM;
    di_p(2,2,:) = param.M_ppM.*(param.alpha_ppM+param.gamma2_ppM);
    di_p(3,2,:) = param.cpII_i_ppM;
    di_p(4,2,:) = param.gamma2_ppM.*param.cpII_i_ppM;
    % third column
    di_p(1,3,:) = 2*param.mu_fr_ppP+param.lambda_ppP+...
        param.alpha_ppP.*param.M_ppP.*param.gamma2_ppP;
    di_p(2,3,:) = param.M_ppP.*(param.alpha_ppP+param.gamma2_ppP);
    di_p(3,3,:) = -param.cpII_i_ppP;
    di_p(4,3,:) = -param.gamma2_ppP.*param.cpII_i_ppP;
    % fourth column
    di_p(1,4,:) = 2*param.mu_fr_ppP+param.lambda_ppP+...
        param.alpha_ppP.*param.M_ppP.*param.gamma1_ppP;
    di_p(2,4,:) = param.M_ppP.*(param.alpha_ppP+param.gamma1_ppP);
    di_p(3,4,:) = -param.cpI_i_ppP;
    di_p(4,4,:) = -param.gamma1_ppP.*param.cpI_i_ppP;

    % Need to invert di to find d. di has large variations of magnitude across 
    % rows but not much within rows. Preconditioning by a diagonal factor to
    % make the maximum  magnitude in each row 1 before inverting greatly 
    % reduces the condition number and in some sense improves the accuracy 
    % of the inverse: in particular, d*di is noticeably closer to the identity,
    % compared to inversion without preconditioning.
   
    % PRECONDITIONING ABANDONED FOR THE TIME BEING. THESE LINES OF CODE
    % ARE FROM THE CONSTANT PHYSICAL PARAMETER CASE
    
    %pre = diag(1.0 ./ sum(abs(di_p), 2));
    %pre = eye(4);
    % fprintf('%.15e\n', diag(pre));
    %d_p = inv(pre*di_p)*pre;

    d_pp = zeros(4,4,num_pp_face_node);
    for fn = 1:num_pp_face_node
        d_pp(:,:,fn) = inv(di_p(:,:,fn));
    end
    
    % save required entries

    param.d11_pp = squeeze(d_pp(1,1,:));
    param.d12_pp = squeeze(d_pp(1,2,:));
    param.d13_pp = squeeze(d_pp(1,3,:));
    param.d14_pp = squeeze(d_pp(1,4,:));
    param.d21_pp = squeeze(d_pp(2,1,:));
    param.d22_pp = squeeze(d_pp(2,2,:));
    param.d23_pp = squeeze(d_pp(2,3,:));
    param.d24_pp = squeeze(d_pp(2,4,:));  

    % maybe for later thought
    % % maximum absolute value of off-diagonal terms of di*d
    % maxdd1 = max(max(abs(di*d - diag(diag(di*d)))));
    % maxdd2 = max(max(abs(d*di - diag(diag(d*di)))));
    % fprintf('max(max(abs(di*d - diag(diag(di*d))))) = %e\n', maxdd1);
    % fprintf('max(max(abs(d*di - diag(diag(d*di))))) = %e\n', maxdd2);

    % r1 eigenvector from \S3.5, 
    % this is used to build flux for PI wave
    param.r1_p_pp = zeros(num_pp_face_node, param.Nfields);
    
    param.r1_p_pp(:, 1) = nx(param.pp_face_node).^2;
    param.r1_p_pp(:, 2) = ny(param.pp_face_node).^2;
    param.r1_p_pp(:, 3) = nz(param.pp_face_node).^2;
    param.r1_p_pp(:, 4) = nx(param.pp_face_node).*ny(param.pp_face_node);
    param.r1_p_pp(:, 5) = ny(param.pp_face_node).*nz(param.pp_face_node);
    param.r1_p_pp(:, 6) = nx(param.pp_face_node).*nz(param.pp_face_node);
    param.r1_p_pp(:, 7) = -param.gamma1_ppM;
    param.r1_p_pp(:, 8) = param.cpI_i_ppM.*nx(param.pp_face_node);
    param.r1_p_pp(:, 9) = param.cpI_i_ppM.*ny(param.pp_face_node);
    param.r1_p_pp(:,10) = param.cpI_i_ppM.*nz(param.pp_face_node);
    param.r1_p_pp(:,11) = param.gamma1_ppM.*param.cpI_i_ppM.*nx(param.pp_face_node);
    param.r1_p_pp(:,12) = param.gamma1_ppM.*param.cpI_i_ppM.*ny(param.pp_face_node);
    param.r1_p_pp(:,13) = param.gamma1_ppM.*param.cpI_i_ppM.*nz(param.pp_face_node);

    % r4 eigenvector from \S3.5, 
    % this is used to build flux for PII wave

    param.r4_p_pp = zeros(num_pp_face_node, param.Nfields);

    param.r4_p_pp(:, 1) = nx(param.pp_face_node).^2;
    param.r4_p_pp(:, 2) = ny(param.pp_face_node).^2;
    param.r4_p_pp(:, 3) = nz(param.pp_face_node).^2;
    param.r4_p_pp(:, 4) = nx(param.pp_face_node).*ny(param.pp_face_node);
    param.r4_p_pp(:, 5) = ny(param.pp_face_node).*nz(param.pp_face_node);
    param.r4_p_pp(:, 6) = nx(param.pp_face_node).*nz(param.pp_face_node);
    param.r4_p_pp(:, 7) = -param.gamma2_ppM;
    param.r4_p_pp(:, 8) = param.cpII_i_ppM.*nx(param.pp_face_node);
    param.r4_p_pp(:, 9) = param.cpII_i_ppM.*ny(param.pp_face_node);
    param.r4_p_pp(:,10) = param.cpII_i_ppM.*nz(param.pp_face_node);
    param.r4_p_pp(:,11) = param.gamma2_ppM.*param.cpII_i_ppM.*nx(param.pp_face_node);
    param.r4_p_pp(:,12) = param.gamma2_ppM.*param.cpII_i_ppM.*ny(param.pp_face_node);
    param.r4_p_pp(:,13) = param.gamma2_ppM.*param.cpII_i_ppM.*nz(param.pp_face_node);

end % if num_pp_face_node > 0

% assemble elastic stiffness and mass matrices on each elastic element. 
% These are Nfields_e square matrices; the third index in the array is the 
% element number within the list of elastic elements.

if num_elastic_elt > 0

    A_e = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);

    A_e( 7, 1, :) = 2*param.mu_ela()+param.lambda_ela;
    A_e( 7, 2, :) = param.lambda_ela;
    A_e( 7, 3, :) = param.lambda_ela;
    A_e( 8, 4, :) = 2*param.mu_ela;
    A_e( 9, 6, :) = 2*param.mu_ela;
    A_e( 1, 7, :) = 1;
    A_e( 4, 8, :) = 0.5;
    A_e( 6, 9, :) = 0.5;
    A_e = -A_e; 
    %param.A_e = A_e;

    B_e = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);
    B_e( 8, 1, :) = param.lambda_ela;
    B_e( 8, 2, :) = 2*param.mu_ela+param.lambda_ela;
    B_e( 8, 3, :) = param.lambda_ela;
    B_e( 7, 4, :) = 2*param.mu_ela;
    B_e( 9, 5, :) = 2*param.mu_ela;
    B_e( 4, 7, :) = 0.5;
    B_e( 2, 8, :) = 1;
    B_e( 5, 9, :) = 0.5;
    B_e = -B_e;
    %param.B_e = B_e;

    C_e = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);
    C_e( 9, 1, :) = param.lambda_ela;
    C_e( 9, 2, :) = param.lambda_ela;
    C_e( 9, 3, :) = 2*param.mu_ela+param.lambda_ela;
    C_e( 8, 5, :) = 2*param.mu_ela;    
    C_e( 7, 6, :) = 2*param.mu_ela;    
    C_e(3,  9, :) = 1;    
    C_e( 5, 8, :) = 0.5;    
    C_e( 6, 7, :) = 0.5;
    C_e = -C_e;
    %param.C_e = C_e;

    Q_e = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);
    Q_e (1, 1, :) = 1;
    Q_e( 2, 2, :) = 1;
    Q_e( 3, 3, :) = 1;
    Q_e( 4, 4, :) = 1;
    Q_e( 5, 5, :) = 1;
    Q_e( 6, 6, :) = 1;
    Q_e( 7, 7, :) = param.rho_ela;
    Q_e( 8, 8, :) = param.rho_ela;
    Q_e( 9, 9, :) = param.rho_ela;
    param.Q_e = Q_e;
      
    % Now we need Q^{-1}, Q^{-1}A, Q^{-1}B, Q^{-1}C for each elastic element.
    % If we're storing the elastic fields in the same arrays as the
    % poroelastic fields, we need to inflate by making rows 7 and 11:13 
    % zero. Do this by allocating zero-initialised storage at the inflated
    % size, computing each Q^{-1}A etc. individually and inserting into the
    % correct rows and columns - this is the role of param.elastic_fields
    
    param.Qi_e  = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    param.QiA_e = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    param.QiB_e = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    param.QiC_e = zeros(param.Nfields, param.Nfields, num_elastic_elt);
        
    for elt = 1:num_elastic_elt
        param.Qi_e (param.elastic_fields, param.elastic_fields, elt) = inv(Q_e(:,:,elt));
        param.QiA_e(param.elastic_fields, param.elastic_fields, elt) = Q_e(:,:,elt) \ A_e(:,:,elt);
        param.QiB_e(param.elastic_fields, param.elastic_fields, elt) = Q_e(:,:,elt) \ B_e(:,:,elt);
        param.QiC_e(param.elastic_fields, param.elastic_fields, elt) = Q_e(:,:,elt) \ C_e(:,:,elt);
    end

end % if num_elastic_elt > 0

if num_ee_face_node > 0

    % flux weights on each elastic-elastic face, refer eqs (106)-(109) in paper. 
    % This builds unnumbered matrix in \S 3.4. It's not given a name there, 
    % but its inverse is d so call it di. It is a 2x2 matrix on each
    % element; the third index is the face node number within
    % elastic-elastic face nodes.
    
    di_ee = zeros(2, 2, num_ee_face_node);
    
    % first column, internal values M
    di_ee(1,1,:) = 2*param.mu_eeM+param.lambda_eeM;
    di_ee(2,1,:) = param.cp_eeM;
    % second column, external values P
    di_ee(1,2,:) = 2*param.mu_eeP+param.lambda_eeP;
    di_ee(2,2,:) = -param.cp_eeP;

    % Need to invert di to find d. Only need the first row, which we
    % compute explicitly using the standard formula.
    % Note that det_di_e is 1 by 1 by num_ee_face_node, as are the
    % di_e(2,2,:), etc., hence use of squeeze
    
    det_di_ee = di_ee(1,1,:).*di_ee(2,2,:) - di_ee(1,2,:).*di_ee(2,1,:);
    param.d11_ee = squeeze( di_ee(2,2,:) ./  det_di_ee);
    param.d12_ee = squeeze(-di_ee(1,2,:) ./  det_di_ee);
%keyboard    
    % ALL THIS PRECONDITIONING CODE LEFT AS-IS IN THE CONSTANT PHYSICAL
    % PARAMETER CASE. REVISIT IF REQUIRED.
    
    % di has large variations of magnitude across 
    % rows but not much within rows. Preconditioning by a diagonal factor to
    % make the maximum  magnitude in each row 1 before inverting greatly 
    % reduces the condition number and in some sense improves the accuracy 
    % of the inverse: in particular, d*di is noticeably closer to the identity,
    % compared to inversion without preconditioning.

    %pre = diag(1.0 ./ sum(abs(di), 2));
    %pre = eye(2);
    % fprintf('%.15e\n', diag(pre));
    %d_e = inv(pre*di_e)*pre;    
    
    % save required entries
    %param.d11_e = d_e(1,1);
    %param.d12_e = d_e(1,2);

    % maybe for later thought
    % % maximum absolute value of off-diagonal terms of di*d
    % maxdd1 = max(max(abs(di*d - diag(diag(di*d)))));
    % maxdd2 = max(max(abs(d*di - diag(diag(d*di)))));
    % fprintf('max(max(abs(di*d - diag(diag(di*d))))) = %e\n', maxdd1);
    % fprintf('max(max(abs(d*di - diag(diag(d*di))))) = %e\n', maxdd2);

    % END PRECONDITIONING CODE FROM CONSTANT PHYSICAL PARAMETER CASE
        
    % r1_ie_ee inflated eigenvector from \S3.5, 
    % this is used to build flux for PI wave
    % if necessary, these are inflated with zeros in the fluid locations 
    % to match the poroelastic case
    param.r1_e_ee = zeros(num_ee_face_node, param.Nfields);
    param.r1_e_ee(:, param.fld_e11) = nx(param.ee_face_node).^2;
    param.r1_e_ee(:, param.fld_e22) = ny(param.ee_face_node).^2;
    param.r1_e_ee(:, param.fld_e33) = nz(param.ee_face_node).^2;
    param.r1_e_ee(:, param.fld_e12) = nx(param.ee_face_node).*ny(param.ee_face_node);
    param.r1_e_ee(:, param.fld_e23) = ny(param.ee_face_node).*nz(param.ee_face_node);
    param.r1_e_ee(:, param.fld_e13) = nx(param.ee_face_node).*nz(param.ee_face_node);
    param.r1_e_ee(:, param.fld_vx)  = param.cp_eeM.*nx(param.ee_face_node);
    param.r1_e_ee(:, param.fld_vy)  = param.cp_eeM.*ny(param.ee_face_node);
    param.r1_e_ee(:, param.fld_vz)  = param.cp_eeM.*nz(param.ee_face_node);

end % if num_ee_face_node > 0
    
% flux weights for interfaces between elastic and poroelastic materials

if num_pe_face_node > 0

    di_pe = zeros(3, 3, num_pe_face_node);
    
    % first column
    di_pe(1,1,:) = 2*param.mu_peM+param.lambda_peM;
    di_pe(2,1,:) = param.cp_peM;
    di_pe(3,1,:) = 0;
    % second column
    di_pe(1,2,:) = 2*param.mu_fr_peP+param.lambda_peP+...
        param.alpha_peP.*param.M_peP.*param.gamma2_peP;
    di_pe(2,2,:) = -param.cpII_i_peP;
    di_pe(3,2,:) = param.gamma2_peP.*param.cpII_i_peP;
    % third column
    di_pe(1,3,:) = 2*param.mu_fr_peP+param.lambda_peP+...
        param.alpha_peP.*param.M_peP.*param.gamma1_peP;
    di_pe(2,3,:) = -param.cpI_i_peP;
    di_pe(3,3,:) = param.gamma1_peP.*param.cpI_i_peP;

    % THIS PRECONDITIONING CODE ABANDONED FOR THE TIME BEING. DATES FROM
    % FIXED PHYSICAL PARAMETER MODEL
    
    %pre = diag(1.0 ./ sum(abs(di), 2));
    %pre = eye(3);
    % fprintf('%.15e\n', diag(pre));
    %d_pe = inv(pre*di_pe)*pre;

    % Now invert each face nodes's di_pe matrix to give a corresponding
    % d_pe matrix
    
    d_pe = zeros(3, 3, num_pe_face_node);
    for fn = 1:num_pe_face_node
        d_pe(:,:,fn) = inv(di_pe(:,:,fn));
    end
    % save required entries
    param.d11_pe = squeeze(d_pe(1,1,:));
    param.d12_pe = squeeze(d_pe(1,2,:));
    param.d13_pe = squeeze(d_pe(1,3,:));
    param.d21_pe = squeeze(d_pe(2,1,:));
    param.d22_pe = squeeze(d_pe(2,2,:));
    param.d23_pe = squeeze(d_pe(2,3,:));
    param.d31_pe = squeeze(d_pe(3,1,:));
    param.d32_pe = squeeze(d_pe(3,2,:));
    param.d33_pe = squeeze(d_pe(3,3,:));

    % r1_ie_pe inflated eigenvector from \S3.5, copied mutatus mutandis 
    % from r1_ie_ee construction for elastic-elastic faces
    % This is elastic so normals are from pe_face_node and parameters
    % are all interior, i.e. M
    
    % Used to build flux for PI wave. Inflated with zeros in the fluid 
    % locations to match the poroelastic case
    param.r1_ie_pe = zeros(num_pe_face_node, param.Nfields);
    
    param.r1_ie_pe(:, 1) = nx(param.pe_face_node).^2;
    param.r1_ie_pe(:, 2) = ny(param.pe_face_node).^2;
    param.r1_ie_pe(:, 3) = nz(param.pe_face_node).^2;
    param.r1_ie_pe(:, 4) = nx(param.pe_face_node).*ny(param.pe_face_node);
    param.r1_ie_pe(:, 5) = ny(param.pe_face_node).*nz(param.pe_face_node);
    param.r1_ie_pe(:, 6) = nx(param.pe_face_node).*nz(param.pe_face_node);
    % row 7 corresponding to zeta not used
    param.r1_ie_pe(:, 8) = param.cp_peM.*nx(param.pe_face_node);
    param.r1_ie_pe(:, 9) = param.cp_peM.*ny(param.pe_face_node);
    param.r1_ie_pe(:,10) = param.cp_peM.*nz(param.pe_face_node);
    % rows 11-13 corresonding to fluid velocity not used

    %***
    
    % r1_p_pe eigenvector from \S3.5, copied mutatus mutandis from r1_p_pp
    % construction for poroelastic-poroelastic faces
    
    % This is poroelastic so normals are ***not*** from pe_face_node 
    % and parameters are all exterior, i.e. P
    % these are constructed using inner pointing normals to poroelastic
    % domain and right propagating waves hence negative signs in velocity 
    % components, consistent with derivation
    
    % Used to build flux for PI wave
    
    param.r1_p_pe = zeros(num_pe_face_node, param.Nfields);
    
    
    param.r1_p_pe(:, 1) = nx(param.pe_face_node).^2;
    param.r1_p_pe(:, 2) = ny(param.pe_face_node).^2;
    param.r1_p_pe(:, 3) = nz(param.pe_face_node).^2;
    param.r1_p_pe(:, 4) = nx(param.pe_face_node).*ny(param.pe_face_node);
    param.r1_p_pe(:, 5) = ny(param.pe_face_node).*nz(param.pe_face_node);
    param.r1_p_pe(:, 6) = nx(param.pe_face_node).*nz(param.pe_face_node);
    param.r1_p_pe(:, 7) = -param.gamma1_peP;
    param.r1_p_pe(:, 8) = -param.cpI_i_peP.*nx(param.pe_face_node);
    param.r1_p_pe(:, 9) = -param.cpI_i_peP.*ny(param.pe_face_node);
    param.r1_p_pe(:,10) = -param.cpI_i_peP.*nz(param.pe_face_node);
    param.r1_p_pe(:,11) = -param.gamma1_peP.*param.cpI_i_peP.*nx(param.pe_face_node);
    param.r1_p_pe(:,12) = -param.gamma1_peP.*param.cpI_i_peP.*ny(param.pe_face_node);
    param.r1_p_pe(:,13) = -param.gamma1_peP.*param.cpI_i_peP.*nz(param.pe_face_node);

    % r4 eigenvector from \S3.5, 
    % this is used to build flux for PII wave

    param.r4_p_pe = zeros(num_pe_face_node, param.Nfields);
    
    param.r4_p_pe(:, 1) = nx(param.pe_face_node).^2;
    param.r4_p_pe(:, 2) = ny(param.pe_face_node).^2;
    param.r4_p_pe(:, 3) = nz(param.pe_face_node).^2;
    param.r4_p_pe(:, 4) = nx(param.pe_face_node).*ny(param.pe_face_node);
    param.r4_p_pe(:, 5) = ny(param.pe_face_node).*nz(param.pe_face_node);
    param.r4_p_pe(:, 6) = nx(param.pe_face_node).*nz(param.pe_face_node);
    param.r4_p_pe(:, 7) = -param.gamma2_peP;
    param.r4_p_pe(:, 8) = -param.cpII_i_peP.*nx(param.pe_face_node);
    param.r4_p_pe(:, 9) = -param.cpII_i_peP.*ny(param.pe_face_node);
    param.r4_p_pe(:,10) = -param.cpII_i_peP.*nz(param.pe_face_node);
    param.r4_p_pe(:,11) = -param.gamma2_peP.*param.cpII_i_peP.*nx(param.pe_face_node);
    param.r4_p_pe(:,12) = -param.gamma2_peP.*param.cpII_i_peP.*ny(param.pe_face_node);
    param.r4_p_pe(:,13) = -param.gamma2_peP.*param.cpII_i_peP.*nz(param.pe_face_node);

% this code is for poroelastic side of elastic/poroelastic interface when
% the poroelatic element is treated as an interior element

    di_ep = zeros(3, 3, num_pe_face_node);
    
    % CHECK SIMON: notationally param.mu_fr_peP shuld be param.mu_fr_epM
    % first column
    di_ep(1,1,:) = 2*param.mu_fr_epM+param.lambda_epM+...
        param.alpha_epM.*param.M_epM.*param.gamma1_epM;
    di_ep(2,1,:) = param.cpI_i_epM;
    di_ep(3,1,:) = param.gamma1_epM.*param.cpI_i_epM;
    % second column
    di_ep(1,2,:) = 2*param.mu_fr_epM+param.lambda_epM+...
        param.alpha_epM.*param.M_epM.*param.gamma2_epM;
    di_ep(2,2,:) = param.cpII_i_epM;
    di_ep(3,2,:) = param.gamma2_epM.*param.cpII_i_epM;
    % third column
    di_ep(1,3,:) = 2*param.mu_epP+param.lambda_epP;
    di_ep(2,3,:) = -param.cp_epP;
    di_ep(3,3,:) = 0;
    
    % Now invert each face nodes's di_ep matrix to give a corresponding
    % d_ep matrix
    
    d_ep = zeros(3, 3, num_pe_face_node);
    for fn = 1:num_pe_face_node
        d_ep(:,:,fn) = inv(di_ep(:,:,fn));
    end
    % save required entries
    param.d11_ep = squeeze(d_ep(1,1,:));
    param.d12_ep = squeeze(d_ep(1,2,:));
    param.d13_ep = squeeze(d_ep(1,3,:));
    param.d21_ep = squeeze(d_ep(2,1,:));
    param.d22_ep = squeeze(d_ep(2,2,:));
    param.d23_ep = squeeze(d_ep(2,3,:));
% UNUSED    
%     param.d31_ep = squeeze(d_ep(3,1,:));
%     param.d32_ep = squeeze(d_ep(3,2,:));
%     param.d33_ep = squeeze(d_ep(3,3,:));

% Used to build flux for PI wave
    
    param.r1_p_ep = zeros(num_pe_face_node, param.Nfields);
    
    % CHECK SIMON: gamma1_peP SHOULD BE gamma1_epM, also
    % nx(param.pe_face_node) should be replaced by  nx(param.ep_face_node)
    % = - nx(param.pe_face_node), so outward pointing relative to poro
    % element
    % SHOULD HAVE SAME STRUCTURE AS r1_p_pp
    
    param.r1_p_ep(:, 1) = nx(param.ep_face_node).^2;
    param.r1_p_ep(:, 2) = ny(param.ep_face_node).^2;
    param.r1_p_ep(:, 3) = nz(param.ep_face_node).^2;
    param.r1_p_ep(:, 4) = nx(param.ep_face_node).*ny(param.ep_face_node);
    param.r1_p_ep(:, 5) = ny(param.ep_face_node).*nz(param.ep_face_node);
    param.r1_p_ep(:, 6) = nx(param.ep_face_node).*nz(param.ep_face_node);
    param.r1_p_ep(:, 7) = -param.gamma1_epM;
    param.r1_p_ep(:, 8) = param.cpI_i_epM.*nx(param.ep_face_node);
    param.r1_p_ep(:, 9) = param.cpI_i_epM.*ny(param.ep_face_node);
    param.r1_p_ep(:,10) = param.cpI_i_epM.*nz(param.ep_face_node);
    param.r1_p_ep(:,11) = param.gamma1_epM.*param.cpI_i_epM.*nx(param.ep_face_node);
    param.r1_p_ep(:,12) = param.gamma1_epM.*param.cpI_i_epM.*ny(param.ep_face_node);
    param.r1_p_ep(:,13) = param.gamma1_epM.*param.cpI_i_epM.*nz(param.ep_face_node);

    % r4 eigenvector from \S3.5, 
    % this is used to build flux for PII wave

    param.r4_p_ep = zeros(num_pe_face_node, param.Nfields);
    
    param.r4_p_ep(:, 1) = nx(param.ep_face_node).^2;
    param.r4_p_ep(:, 2) = ny(param.ep_face_node).^2;
    param.r4_p_ep(:, 3) = nz(param.ep_face_node).^2;
    param.r4_p_ep(:, 4) = nx(param.ep_face_node).*ny(param.ep_face_node);
    param.r4_p_ep(:, 5) = ny(param.ep_face_node).*nz(param.ep_face_node);
    param.r4_p_ep(:, 6) = nx(param.ep_face_node).*nz(param.ep_face_node);
    param.r4_p_ep(:, 7) = -param.gamma2_epM;
    param.r4_p_ep(:, 8) = param.cpII_i_epM.*nx(param.ep_face_node);
    param.r4_p_ep(:, 9) = param.cpII_i_epM.*ny(param.ep_face_node);
    param.r4_p_ep(:,10) = param.cpII_i_epM.*nz(param.ep_face_node);
    param.r4_p_ep(:,11) = param.gamma2_epM.*param.cpII_i_epM.*nx(param.ep_face_node);
    param.r4_p_ep(:,12) = param.gamma2_epM.*param.cpII_i_epM.*ny(param.ep_face_node);
    param.r4_p_ep(:,13) = param.gamma2_epM.*param.cpII_i_epM.*nz(param.ep_face_node);
    
end % if num_interface_elt > 0

% END ASSEMBLY OF ELASTIC AND POROELASTIC PARAMETERS

% BEGIN PLANE WAVE SETUP FOR CONVERGENCE CODE
% CONSIDER MOVING TO SEPARATE FILE
% BODGED BY USING THE FIRST ELASTIC OR FIRST POROELASTIC ELEMENT TO COMPUTE
% EIGENVALUES AND EIGENVECTORS

% Used to build analytic plane waves propagating through (1,1,1)
% To do: move this stuff to build_param_plane_wave since only used in
% convergence code. HOWEVER: the computed wavespeeds have generic use,
% write get_low_frequency_wavespeed and get_high_frequency_wavespeed using
% formulae in appendix of paper
% nnx = 1/sqrt(3);
% nny = 1/sqrt(3);
% nnz = 1/sqrt(3); % THIS OVERWRITES MATLAB nnz FUNCTION!

%direction cosines
nnx = param.pwx/sqrt(param.pwx^2+param.pwy^2+param.pwz^2); 
nny = param.pwy/sqrt(param.pwx^2+param.pwy^2+param.pwz^2); 
nnz = param.pwz/sqrt(param.pwx^2+param.pwy^2+param.pwz^2); % THIS OVERWRITES MATLAB nnz FUNCTION!

if num_elastic_elt > 0

    % used to build elastic plane wave in plane_wave_elastic
    Pi = nnx*A_e(:,:,1) + nny*B_e(:,:,1) + nnz*C_e(:,:,1);
    [e_vec,e_val] = eig(Pi, Q_e(:,:,1));
    e_val = diag(e_val);
    % sort eigenvalues and eigenvectors
    [~, ind] = sort(e_val);
    param.e_vec_plane_wave_ela = e_vec(:,ind);
    param.c_plane_wave_ela = e_val(ind);

    % Now inflate if necessary by adding zero rows in positions 7, 10:13

    param.e_vec_plane_wave_e = zeros(param.Nfields, param.Nfields_e);
    param.e_vec_plane_wave_e(param.elastic_fields,:) = param.e_vec_plane_wave_ela;

end % if num_elastic_elt > 0

if num_inviscid_elt > 0

    % used to build inviscid poroelatic plane wave in plane_wave_inviscid 
    Pi = nnx*param.A_p(:,:,1) + nny*param.B_p(:,:,1) + nnz*param.C_p(:,:,1);
    [e_vec,e_val] = eig(Pi, param.Q_p(:,:,1));
    e_val = diag(e_val);
    % sort eigenvalues and eigenvectors
    [~, ind] = sort(e_val);
    param.e_vec_plane_wave_in = e_vec(:,ind);
    param.c_plane_wave_in = e_val(ind);

end % if num_inviscid_elt > 0

if num_lf_elt > 0

    % used to build diffusive low frequency plane wave
    Pi = nnx*param.A_p(:,:,1) + nny*param.B_p(:,:,1) + nnz*param.C_p(:,:,1);
    [e_vec_lf,e_val_lf] = eig(Pi, param.Q_p(:,:,1)-1i/param.omega_sou*D_p(:,:,1));
    e_val_lf = diag(e_val_lf);
    % sort eigenvalues and eigenvectors
    [~, ind] = sort(real(e_val_lf));
    param.e_vec_plane_wave_lf = e_vec_lf(:,ind);
    param.c_plane_wave_lf = e_val_lf(ind);
    param.ph_vel_lf = 1./real(1./e_val_lf(ind));
    
    % for operator-splitting and IMEX
    param.rho_f_lf = param.rho_f(param.lf_elt);
    param.rho_a_lf = param.rho_a(param.lf_elt);
    param.m_lf = param.m(param.lf_elt); % needed for IMEX
    param.stiff_eigenval =  -param.rho_a_lf.*param.eta_lf ...
        ./((param.m_lf.*param.rho_a_lf-param.rho_f_lf.^2).*param.k_lf);
    
end % if num_lf_elt > 0

if num_hf_elt > 0
    % used to build high frequency plane wave
    Pi_hf = nnx*param.A_hf(:,:,1) + nny*param.B_hf(:,:,1) + nnz*param.C_hf(:,:,1);
    [e_vec_hf,e_val_hf] = eig(Pi_hf, param.Q_hf(:,:,1)-1i/param.omega_sou*param.D_hf(:,:,1));
    e_val_hf = diag(e_val_hf);
    % sort eigenvalues and eigenvectors
    [~, ind] = sort(real(e_val_hf));
    param.e_vec_plane_wave_hf = e_vec_hf(:,ind);
    param.c_plane_wave_hf = e_val_hf(ind);
    param.ph_vel_hf = 1./real(1./e_val_hf(ind));
end % if num_hf_elt > 0

% END PLANE WAVE SETUP FOR CONVERGENCE CODE
    
return

