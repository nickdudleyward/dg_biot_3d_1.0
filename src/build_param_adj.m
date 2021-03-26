function param = build_param_adj(param)
% Builds data structure with material parameters for each element for
% elastic and poroelastic cases
% this is only implemented for inviscid case

Globals3D;

num_elastic_elt     = numel(param.elastic_elt);
num_poroelastic_elt = numel(param.poroelastic_elt);
num_hf_elt          = numel(param.hf_elt);
num_lf_elt          = numel(param.lf_elt);
num_inviscid_elt    = numel(param.inviscid_elt);

num_ee_face_node    = numel(param.ee_face_node);
num_pp_face_node    = numel(param.pp_face_node);
num_pe_face_node    = numel(param.pe_face_node);

% assemble poroelastic elastic stiffness and mass matrices on each 
% poroelastic element. These are Nfields square matrices; the third index 
% in the array is the element number within the list of poroelastic 
% elements.

if num_poroelastic_elt > 0

    A_p_adj = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);
    A_p_adj(  8, 1, :) = 1;
    A_p_adj(  9, 4, :) = 1;
    A_p_adj( 10, 6, :) = 1;
    A_p_adj( 11, 7, :) = -1;

    A_p_adj( 1,  8, :) = 2*param.mu_fr+param.lambda;
    A_p_adj( 1, 11, :) = param.M.*param.alpha;
    A_p_adj( 2,  8, :) = param.lambda;
    A_p_adj( 2, 11, :) = param.M.*param.alpha;
    A_p_adj( 3,  8, :) = param.lambda;
    A_p_adj( 3, 11, :) = param.M.*param.alpha;
    A_p_adj( 4,  9, :) = param.mu_fr;
    A_p_adj( 6, 10, :) = param.mu_fr;
    A_p_adj( 7,  8, :) = -param.alpha.*param.M;
    A_p_adj( 7, 11, :) = -param.M;
    
    A_p_adj = -A_p_adj; % this is to make consistent with the hyp system formulation d/dt + d/dx(aq) +... =  sources
    param.A_p_adj = A_p_adj;
    
    B_p_adj = zeros(param.Nfields,param.Nfields, num_poroelastic_elt);
    B_p_adj( 8, 4, :) = 1;
    B_p_adj( 9, 2, :) = 1;
    B_p_adj(10, 5, :) = 1;
    B_p_adj(12, 7, :) = -1;

    
    B_p_adj( 1, 9, :) = param.lambda;
    B_p_adj( 1,12, :) = param.M.*param.alpha;
    B_p_adj( 2, 9, :) = 2*param.mu_fr+param.lambda;
    B_p_adj( 2,12, :) = param.M.*param.alpha;
    B_p_adj( 3, 9, :) = param.lambda;
    B_p_adj( 3,12, :) = param.M.*param.alpha;
    B_p_adj( 4, 8, :) = param.mu_fr;
    B_p_adj( 5,10, :) = param.mu_fr;
    B_p_adj( 7, 9, :) = -param.alpha.*param.M;
    B_p_adj( 7,12, :) = -param.M;
    
    B_p_adj = -B_p_adj;
    param.B_p_adj = B_p_adj;

    C_p_adj = zeros(param.Nfields,param.Nfields, num_poroelastic_elt);
    C_p_adj( 8, 6, :) = 1;
    C_p_adj( 9, 5, :) = 1; 
    C_p_adj(10, 3, :) = 1;    
    C_p_adj(13, 7, :) = -1;
    
    C_p_adj( 1, 10, :) = param.lambda;
    C_p_adj( 1, 13, :) = param.M.*param.alpha;    
    C_p_adj( 2, 10, :) = param.lambda;
    C_p_adj( 2, 13, :) = param.M.*param.alpha;    
    C_p_adj( 3, 10, :) = 2*param.mu_fr+param.lambda;
    C_p_adj( 3, 13, :) = param.M.*param.alpha;    
    C_p_adj( 5,  9, :) = param.mu_fr;    
    C_p_adj( 6,  8, :) = param.mu_fr;   
    C_p_adj( 7, 10, :) = -param.alpha.*param.M;
    C_p_adj( 7, 13, :) = -param.M;
    
    C_p_adj = -C_p_adj;
    param.C_p_adj = C_p_adj;

%     % Diffusivity for low frequency case
% 
%     D_p = zeros(param.Nfields, param.Nfields, num_lf_elt);
%     D_p(11,11, :) = -param.eta_lf./param.k_lf;
%     D_p(12,12, :) = -param.eta_lf./param.k_lf;
%     D_p(13,13, :) = -param.eta_lf./param.k_lf;
%     param.D_p = D_p;
    
    % poroelastic mass matrix 
    % this is self-adjoint, for now we will keep it
    Q_p_adj = param.Q_p;

    % Now we need Q^{-1}A, Q^{-1}B, Q^{-1}C for each poroelastic element. 
    % Handle each element individually.
    
    % Also need Q^{-1} itself, for the LF dissipative case and for any
    % poroelastic (not necessarily LF) element containing a source.
    % This is a pain, and is done very inefficiently here. We compute Qi_p
    % for all poroelastic elements and for each LF element take a copy.
    % This still makes finding the Qi_p for the source element a pain
    % because param.elm_sou is an index into all elements but param.Qi_p
    % contains only poroelastic elements

    param.Qi_p_adj   = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);    
    param.QiA_p_adj  = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);
    param.QiB_p_adj  = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);
    param.QiC_p_adj  = zeros(param.Nfields, param.Nfields, num_poroelastic_elt);

    param.Qi_lf_adj  = zeros(param.Nfields, param.Nfields, num_lf_elt);

    elt_lf = 1;
    for elt = 1:num_poroelastic_elt
        param.Qi_p_adj(:,:,elt) = inv(Q_p_adj(:,:,elt));
        param.QiA_p_adj(:,:,elt) =  Q_p_adj(:,:,elt) \ A_p_adj(:,:,elt);
        param.QiB_p_adj(:,:,elt) =  Q_p_adj(:,:,elt) \ B_p_adj(:,:,elt);
        param.QiC_p_adj(:,:,elt) =  Q_p_adj(:,:,elt) \ C_p_adj(:,:,elt);
        % if we have a LF eleemnt, copy Qi to Qi_lf
        if param.elt_to_model_type(elt) == param.LOW_FREQUENCY
            param.Qi_lf_adj(:,:,elt_lf) = param.Qi_p_adj(:,:,elt);
            elt_lf  = elt_lf + 1;
        end
    end % for elt = 1:num_poroelastic_elt
    
end % if num_poroelastic_elt > 0

% sanity checks

% fprintf('if assembled correctly should return zero matrix')
% temp = param.A_p(:,:,1);
% temp(4:6,:) = 2*temp(4:6,:);
% temp(9:10,:) = temp(9:10,:)/2;
% temp'-param.A_p_adj(:,:,1)
% 
% temp = param.B_p(:,:,1);
% temp(4:5,:) = 2*temp(4:5,:);
% temp([8 10],:) = temp([8 10],:)/2;
% temp'-param.B_p_adj(:,:,1)
% 
% temp = param.C_p(:,:,1);
% temp(5:6,:) = 2*temp(5:6,:);
% temp(8:9,:) = temp(8:9,:)/2;
% temp'-param.C_p_adj(:,:,1)
% 
% fprintf('should print out same eigenvalues')
% eig(param.QiA_p(:,:,1))
% eig(param.QiA_p_adj(:,:,1))
% eig(param.QiB_p(:,:,1))
% eig(param.QiB_p_adj(:,:,1))
% eig(param.QiC_p(:,:,1))
% eig(param.QiC_p_adj(:,:,1))

%--------------------------------------------------------------------------
% HIGH FREQUENCY TO DO

if num_hf_elt > 0
    
    fprintf('To be coded')
% 
%     % Storage for high-frequency versions, with three extra rows and
%     % columns to account for memory variables.
%     
%     % These are used for analytic plane waves and QiD is used in
%     % dissipation_hf
%     
%     param.A_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
%     param.B_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
%     param.C_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
%     param.D_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
%     param.Q_hf   = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
%     param.QiD_hf = zeros(param.Nfields+3, param.Nfields+3, num_hf_elt);
%     
%     % Extend A_p, B_p, C_p to high-frequency versions adding three zero 
%     % rows and columns and the bottom and right   
%     
%     param.A_hf(1:param.Nfields, 1:param.Nfields, :) = A_p(:,:,param.hf_within_poro);
%     param.B_hf(1:param.Nfields, 1:param.Nfields, :) = B_p(:,:,param.hf_within_poro);
%     param.C_hf(1:param.Nfields, 1:param.Nfields, :) = C_p(:,:,param.hf_within_poro);
% 
%     % high frequency mass matrix, extended to include memory variables
%     param.Q_hf(1:param.Nfields, 1:param.Nfields, :) = Q_p(:,:,param.hf_within_poro);
%     param.Q_hf(14,14,:) = 1;
%     param.Q_hf(15,15,:) = 1;
%     param.Q_hf(16,16,:) = 1;
%     param.Q_hf(14,11,:) = -(param.taue_hf./param.taus_hf-1);
%     param.Q_hf(15,12,:) = -(param.taue_hf./param.taus_hf-1);
%     param.Q_hf(16,13,:) = -(param.taue_hf./param.taus_hf-1);
% 
%     % high frequency diffusivity matrix , extended to include memory variables
%     param.D_hf = zeros(param.Nfields, param.Nfields, num_hf_elt);
%     param.D_hf(11,11,:)   = -param.eta_hf./param.k_hf;
%     param.D_hf(12,12,:)   = -param.eta_hf./param.k_hf;
%     param.D_hf(13,13,:)   = -param.eta_hf./param.k_hf;
%     param.D_hf(11,14,:)   = -param.eta_hf./param.k_hf;
%     param.D_hf(12,15,:)   = -param.eta_hf./param.k_hf;
%     param.D_hf(13,16,:)   = -param.eta_hf./param.k_hf;
%     param.D_hf(14,14,:)   = -1/param.taus_hf;
%     param.D_hf(15,15,:)   = -1/param.taus_hf;
%     param.D_hf(16,16,:)   = -1/param.taus_hf;
% 
%     for elt=1:num_hf_elt
%         param.QiD_hf(:,:,elt) = param.Q_hf(:,:,elt) \ param.D_hf(:,:,elt);
%     end
%     
end % if num_hf_elt > 0

%--------------------------------------------------------------------------

if num_pp_face_node > 0
    
    %flux stuff, refer eqs (106)-(109) in paper. 
    
    % first note that the d weights have already been defined and are the
    % identical for the adjoint fluxes
    
    
    % param.d11_pp = squeeze(d_pp(1,1,:));
    % param.d12_pp = squeeze(d_pp(1,2,:));
    % param.d13_pp = squeeze(d_pp(1,3,:));
    % param.d14_pp = squeeze(d_pp(1,4,:));
    % param.d21_pp = squeeze(d_pp(2,1,:));
    % param.d22_pp = squeeze(d_pp(2,2,:));
    % param.d23_pp = squeeze(d_pp(2,3,:));
    % param.d24_pp = squeeze(d_pp(2,4,:));  
    
    

    % r1 eigenvector from \S3.5, 
    % this is used to build flux for PI wave
    param.r1_p_pp_adj = zeros(num_pp_face_node, param.Nfields);
    
    param.r1_p_pp_adj(:, 1) = param.lambda_ppM + 2*param.mu_fr_ppM.*nx(param.pp_face_node).^2  + param.alpha_ppM.*param.gamma1_ppM.*param.M_ppM;
    param.r1_p_pp_adj(:, 2) = param.lambda_ppM + 2*param.mu_fr_ppM.*ny(param.pp_face_node).^2  + param.alpha_ppM.*param.gamma1_ppM.*param.M_ppM;
    param.r1_p_pp_adj(:, 3) = param.lambda_ppM + 2*param.mu_fr_ppM.*nz(param.pp_face_node).^2  + param.alpha_ppM.*param.gamma1_ppM.*param.M_ppM;  
    
    param.r1_p_pp_adj(:, 4) = 2*param.mu_fr_ppM.*nx(param.pp_face_node).*ny(param.pp_face_node);
    param.r1_p_pp_adj(:, 5) = 2*param.mu_fr_ppM.*ny(param.pp_face_node).*nz(param.pp_face_node);
    param.r1_p_pp_adj(:, 6) = 2*param.mu_fr_ppM.*nx(param.pp_face_node).*nz(param.pp_face_node);
    
    
    param.r1_p_pp_adj(:, 7) = -param.alpha_ppM.*param.M_ppM - param.gamma1_ppM.*param.M_ppM;
    
    param.r1_p_pp_adj(:, 8) = param.cpI_i_ppM.*nx(param.pp_face_node);
    param.r1_p_pp_adj(:, 9) = param.cpI_i_ppM.*ny(param.pp_face_node);
    param.r1_p_pp_adj(:,10) = param.cpI_i_ppM.*nz(param.pp_face_node);
    param.r1_p_pp_adj(:,11) = param.gamma1_ppM.*param.cpI_i_ppM.*nx(param.pp_face_node);
    param.r1_p_pp_adj(:,12) = param.gamma1_ppM.*param.cpI_i_ppM.*ny(param.pp_face_node);
    param.r1_p_pp_adj(:,13) = param.gamma1_ppM.*param.cpI_i_ppM.*nz(param.pp_face_node);

    % r4 eigenvector from \S3.5, 
    % this is used to build flux for PII wave

    param.r4_p_pp_adj = zeros(num_pp_face_node, param.Nfields);

    param.r4_p_pp_adj(:, 1) = param.lambda_ppM + 2*param.mu_fr_ppM.*nx(param.pp_face_node).^2  + param.alpha_ppM.*param.gamma2_ppM.*param.M_ppM;
    param.r4_p_pp_adj(:, 2) = param.lambda_ppM + 2*param.mu_fr_ppM.*ny(param.pp_face_node).^2  + param.alpha_ppM.*param.gamma2_ppM.*param.M_ppM;
    param.r4_p_pp_adj(:, 3) = param.lambda_ppM + 2*param.mu_fr_ppM.*nz(param.pp_face_node).^2  + param.alpha_ppM.*param.gamma2_ppM.*param.M_ppM;
    
    param.r4_p_pp_adj(:, 4) = 2*param.mu_fr_ppM.*nx(param.pp_face_node).*ny(param.pp_face_node);
    param.r4_p_pp_adj(:, 5) = 2*param.mu_fr_ppM.*ny(param.pp_face_node).*nz(param.pp_face_node);
    param.r4_p_pp_adj(:, 6) = 2*param.mu_fr_ppM.*nx(param.pp_face_node).*nz(param.pp_face_node);
    
    param.r4_p_pp_adj(:, 7) = -param.alpha_ppM.*param.M_ppM -param.gamma2_ppM.*param.M_ppM;
    
    param.r4_p_pp_adj(:, 8) = param.cpII_i_ppM.*nx(param.pp_face_node);
    param.r4_p_pp_adj(:, 9) = param.cpII_i_ppM.*ny(param.pp_face_node);
    param.r4_p_pp_adj(:,10) = param.cpII_i_ppM.*nz(param.pp_face_node);
    param.r4_p_pp_adj(:,11) = param.gamma2_ppM.*param.cpII_i_ppM.*nx(param.pp_face_node);
    param.r4_p_pp_adj(:,12) = param.gamma2_ppM.*param.cpII_i_ppM.*ny(param.pp_face_node);
    param.r4_p_pp_adj(:,13) = param.gamma2_ppM.*param.cpII_i_ppM.*nz(param.pp_face_node);

end % if num_pp_face_node > 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assemble elastic stiffness and mass matrices on each elastic element. 
% These are Nfields_e square matrices; the third index in the array is the 
% element number within the list of elastic elements.

if num_elastic_elt > 0

    A_e_adj = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);

    A_e_adj( 1, 7, :) = 2*param.mu_ela+param.lambda_ela;
    A_e_adj( 2, 7, :) = param.lambda_ela;
    A_e_adj( 3, 7, :) = param.lambda_ela;
    A_e_adj( 4, 8, :) = param.mu_ela;
    A_e_adj( 6, 9, :) = param.mu_ela;
    A_e_adj( 7, 1, :) = 1;
    A_e_adj( 8, 4, :) = 1;
    A_e_adj( 9, 6, :) = 1;
    A_e_adj = -A_e_adj; 
    %param.A_e_adj = A_e_adj;

    B_e_adj = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);
    B_e_adj( 1, 8, :) = param.lambda_ela;
    B_e_adj( 2, 8, :) = 2*param.mu_ela+param.lambda_ela;
    B_e_adj( 3, 8, :) = param.lambda_ela;
    B_e_adj( 4, 7, :) = param.mu_ela;
    B_e_adj( 5, 9, :) = param.mu_ela;
    B_e_adj( 7, 4, :) = 1;
    B_e_adj( 8, 2, :) = 1;
    B_e_adj( 9, 5, :) = 1;
    B_e_adj = -B_e_adj;
    %param.B_e_adj = B_e_adj;

    C_e_adj = zeros(param.Nfields_e,param.Nfields_e, num_elastic_elt);
    C_e_adj( 1, 9, :) = param.lambda_ela;
    C_e_adj( 2, 9, :) = param.lambda_ela;
    C_e_adj( 3, 9, :) = 2*param.mu_ela+param.lambda_ela;
    C_e_adj( 5, 8, :) = param.mu_ela;    
    C_e_adj( 6, 7, :) = param.mu_ela;    
    C_e_adj( 7, 6, :) = 1;    
    C_e_adj( 8, 5, :) = 1;    
    C_e_adj( 9, 3, :) = 1;
    C_e_adj = -C_e_adj;
    %param.C_e_adj = C_e_adj;
    
    % poroelastic mass matrix 
    % this is self-adjoint, for now we will keep it
    Q_e_adj = param.Q_e;

    % Now we need Q^{-1}, Q^{-1}A, Q^{-1}B, Q^{-1}C for each elastic element.
    % If we're storing the elastic fields in the same arrays as the
    % poroelastic fields, we need to inflate by making rows 7 and 11:13 
    % zero. Do this by allocating zero-initialised storage at the inflated
    % size, computing each Q^{-1}A etc. individually and inserting into the
    % correct rows and columns - this is the role of param.elastic_fields
    
    param.Qi_e_adj  = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    param.QiA_e_adj = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    param.QiB_e_adj = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    param.QiC_e_adj = zeros(param.Nfields, param.Nfields, num_elastic_elt);
    for elt = 1:num_elastic_elt
        param.Qi_e_adj (param.elastic_fields, param.elastic_fields, elt) = inv(Q_e_adj(:,:,elt));
        param.QiA_e_adj(param.elastic_fields, param.elastic_fields, elt) = Q_e_adj(:,:,elt) \ A_e_adj(:,:,elt);
        param.QiB_e_adj(param.elastic_fields, param.elastic_fields, elt) = Q_e_adj(:,:,elt) \ B_e_adj(:,:,elt);
        param.QiC_e_adj(param.elastic_fields, param.elastic_fields, elt) = Q_e_adj(:,:,elt) \ C_e_adj(:,:,elt);
    end

end % if num_elastic_elt > 0

% TEST CODE
% fprintf('should print out same eigenvalues')
% eig(param.QiA_ie(:,:,1))
% eig(param.QiA_ie_adj(:,:,1))
% eig(param.QiB_ie(:,:,1))
% eig(param.QiB_ie_adj(:,:,1))
% eig(param.QiC_ie(:,:,1))
% eig(param.QiC_ie_adj(:,:,1))

if num_ee_face_node > 0
        
    % r1_ie_ee inflated eigenvector from \S3.5, 
    % this is used to build flux for PI wave
    % these are inflated with zeros in the fluid locations to match the
    % poroelastic case
    param.r1_e_ee_adj = zeros(num_ee_face_node, param.Nfields);
    param.r1_e_ee_adj(:, param.fld_s11) = param.lambda_eeM + 2*param.mu_eeM.*nx(param.ee_face_node).^2;
    param.r1_e_ee_adj(:, param.fld_s22) = param.lambda_eeM + 2*param.mu_eeM.*ny(param.ee_face_node).^2;
    param.r1_e_ee_adj(:, param.fld_s33) = param.lambda_eeM + 2*param.mu_eeM.*nz(param.ee_face_node).^2;    
    param.r1_e_ee_adj(:, param.fld_s12) = 2*param.mu_eeM.*nx(param.ee_face_node).*ny(param.ee_face_node);
    param.r1_e_ee_adj(:, param.fld_s23) = 2*param.mu_eeM.*ny(param.ee_face_node).*nz(param.ee_face_node);
    param.r1_e_ee_adj(:, param.fld_s13) = 2*param.mu_eeM.*nx(param.ee_face_node).*nz(param.ee_face_node);
    param.r1_e_ee_adj(:, param.fld_vx) = param.cp_eeM.*nx(param.ee_face_node);
    param.r1_e_ee_adj(:, param.fld_vy) = param.cp_eeM.*ny(param.ee_face_node);
    param.r1_e_ee_adj(:, param.fld_vz) = param.cp_eeM.*nz(param.ee_face_node);
    
end % if num_ee_face_node > 0

% flux weights for interfaces between elastic and poroelastic materials

if num_pe_face_node > 0
    
    % r1_ie_pe inflated eigenvector from \S3.5, copied mutatus mutandis
    % from r1_ie_ee construction for elastic-elastic faces
    % This is elastic so normals are from pe_face_node and parameters
    % are all interior, i.e. M
    
    % Used to build flux for PI wave. Inflated with zeros in the fluid
    % locations to match the poroelastic case
    
    param.r1_ie_pe_adj = zeros(num_pe_face_node, param.Nfields);
    
    param.r1_ie_pe_adj(:, 1) = param.lambda_peM + 2*param.mu_peM.*nx(param.pe_face_node).^2;
    param.r1_ie_pe_adj(:, 2) = param.lambda_peM + 2*param.mu_peM.*ny(param.pe_face_node).^2;
    param.r1_ie_pe_adj(:, 3) = param.lambda_peM + 2*param.mu_peM.*nz(param.pe_face_node).^2;
    
    param.r1_ie_pe_adj(:, 4) = 2*param.mu_peM.*nx(param.pe_face_node).*ny(param.pe_face_node);
    param.r1_ie_pe_adj(:, 5) = 2*param.mu_peM.*ny(param.pe_face_node).*nz(param.pe_face_node);
    param.r1_ie_pe_adj(:, 6) = 2*param.mu_peM.*nx(param.pe_face_node).*nz(param.pe_face_node);
    % row 7 corresponding to zeta not used
    param.r1_ie_pe_adj(:, 8) = param.cp_peM.*nx(param.pe_face_node);
    param.r1_ie_pe_adj(:, 9) = param.cp_peM.*ny(param.pe_face_node);
    param.r1_ie_pe_adj(:,10) = param.cp_peM.*nz(param.pe_face_node);
    % rows 11-13 corresonding to fluid velocity not used
    
    % r1_p_pe eigenvector from \S3.5, copied mutatus mutandis from r1_p_pp
    % construction for poroelastic-poroelastic faces
    
    % This is poroelastic so normals are ***not*** from pe_face_node
    % and parameters are all exterior, i.e. P
    % these are constructed using inner pointing normals to poroelastic
    % domain and right propagating waves hence negative signs in velocity
    % components, consistent with derivation
    
    % Used to build flux for PI wave
    
    param.r1_p_pe_adj = zeros(num_pe_face_node, param.Nfields);
    
    param.r1_p_pe_adj(:, 1) = param.lambda_peP + 2*param.mu_fr_peP.*nx(param.pe_face_node).^2  + param.alpha_peP.*param.gamma1_peP.*param.M_peP;
    param.r1_p_pe_adj(:, 2) = param.lambda_peP + 2*param.mu_fr_peP.*ny(param.pe_face_node).^2  + param.alpha_peP.*param.gamma1_peP.*param.M_peP;
    param.r1_p_pe_adj(:, 3) = param.lambda_peP + 2*param.mu_fr_peP.*nz(param.pe_face_node).^2  + param.alpha_peP.*param.gamma1_peP.*param.M_peP;
    
    param.r1_p_pe_adj(:, 4) = 2*param.mu_fr_peP.*nx(param.pe_face_node).*ny(param.pe_face_node);
    param.r1_p_pe_adj(:, 5) = 2*param.mu_fr_peP.*ny(param.pe_face_node).*nz(param.pe_face_node);
    param.r1_p_pe_adj(:, 6) = 2*param.mu_fr_peP.*nx(param.pe_face_node).*nz(param.pe_face_node);
    
    param.r1_p_pe_adj(:, 7) = -param.alpha_peP.*param.M_peP - param.gamma1_peP.*param.M_peP;
    
    param.r1_p_pe_adj(:, 8) = -param.cpI_i_peP.*nx(param.pe_face_node);
    param.r1_p_pe_adj(:, 9) = -param.cpI_i_peP.*ny(param.pe_face_node);
    param.r1_p_pe_adj(:,10) = -param.cpI_i_peP.*nz(param.pe_face_node);
    param.r1_p_pe_adj(:,11) = -param.gamma1_peP.*param.cpI_i_peP.*nx(param.pe_face_node);
    param.r1_p_pe_adj(:,12) = -param.gamma1_peP.*param.cpI_i_peP.*ny(param.pe_face_node);
    param.r1_p_pe_adj(:,13) = -param.gamma1_peP.*param.cpI_i_peP.*nz(param.pe_face_node);
    
    % r4 eigenvector from \S3.5,
    % this is used to build flux for PII wave
    
    param.r4_p_pe_adj = zeros(num_pe_face_node, param.Nfields);
    
    param.r4_p_pe_adj(:, 1) = param.lambda_peP + 2*param.mu_fr_peP.*nx(param.pe_face_node).^2  + param.alpha_peP.*param.gamma2_peP.*param.M_peP;
    param.r4_p_pe_adj(:, 2) = param.lambda_peP + 2*param.mu_fr_peP.*ny(param.pe_face_node).^2  + param.alpha_peP.*param.gamma2_peP.*param.M_peP;
    param.r4_p_pe_adj(:, 3) = param.lambda_peP + 2*param.mu_fr_peP.*nz(param.pe_face_node).^2  + param.alpha_peP.*param.gamma2_peP.*param.M_peP;
    
    param.r4_p_pe_adj(:, 4) = 2*param.mu_fr_peP.*nx(param.pe_face_node).*ny(param.pe_face_node);
    param.r4_p_pe_adj(:, 5) = 2*param.mu_fr_peP.*ny(param.pe_face_node).*nz(param.pe_face_node);
    param.r4_p_pe_adj(:, 6) = 2*param.mu_fr_peP.*nx(param.pe_face_node).*nz(param.pe_face_node);
    
    param.r4_p_pe_adj(:, 7) = -param.alpha_peP.*param.M_peP - param.gamma2_peP.*param.M_peP;
    
    param.r4_p_pe_adj(:, 8) = -param.cpII_i_peP.*nx(param.pe_face_node);
    param.r4_p_pe_adj(:, 9) = -param.cpII_i_peP.*ny(param.pe_face_node);
    param.r4_p_pe_adj(:,10) = -param.cpII_i_peP.*nz(param.pe_face_node);
    param.r4_p_pe_adj(:,11) = -param.gamma2_peP.*param.cpII_i_peP.*nx(param.pe_face_node);
    param.r4_p_pe_adj(:,12) = -param.gamma2_peP.*param.cpII_i_peP.*ny(param.pe_face_node);
    param.r4_p_pe_adj(:,13) = -param.gamma2_peP.*param.cpII_i_peP.*nz(param.pe_face_node);
    
    % this code is for poroelastic side of elastic/poroelastic interface when
    % the poroelatic element is treated as an interior element
    % Used to build flux for PI wave
       
    param.r1_p_ep_adj = zeros(num_pe_face_node, param.Nfields);
    
    % SHOULD HAVE SAME STRUCTURE AS r1_p_pp_adj
    
    param.r1_p_ep_adj(:, 1) = param.lambda_epM + 2*param.mu_fr_epM.*nx(param.ep_face_node).^2  + param.alpha_epM.*param.gamma1_epM.*param.M_epM;
    param.r1_p_ep_adj(:, 2) = param.lambda_epM + 2*param.mu_fr_epM.*ny(param.ep_face_node).^2  + param.alpha_epM.*param.gamma1_epM.*param.M_epM;
    param.r1_p_ep_adj(:, 3) = param.lambda_epM + 2*param.mu_fr_epM.*nz(param.ep_face_node).^2  + param.alpha_epM.*param.gamma1_epM.*param.M_epM;
    
    param.r1_p_ep_adj(:, 4) = 2*param.mu_fr_epM.*nx(param.ep_face_node).*ny(param.ep_face_node);
    param.r1_p_ep_adj(:, 5) = 2*param.mu_fr_epM.*ny(param.ep_face_node).*nz(param.ep_face_node);
    param.r1_p_ep_adj(:, 6) = 2*param.mu_fr_epM.*nx(param.ep_face_node).*nz(param.ep_face_node);
    
    param.r1_p_ep_adj(:, 7) = -param.alpha_epM.*param.M_epM - param.gamma1_epM.*param.M_epM;
    
    param.r1_p_ep_adj(:, 8) = param.cpI_i_epM.*nx(param.ep_face_node);
    param.r1_p_ep_adj(:, 9) = param.cpI_i_epM.*ny(param.ep_face_node);
    param.r1_p_ep_adj(:,10) = param.cpI_i_epM.*nz(param.ep_face_node);
    param.r1_p_ep_adj(:,11) = param.gamma1_epM.*param.cpI_i_epM.*nx(param.ep_face_node);
    param.r1_p_ep_adj(:,12) = param.gamma1_epM.*param.cpI_i_epM.*ny(param.ep_face_node);
    param.r1_p_ep_adj(:,13) = param.gamma1_epM.*param.cpI_i_epM.*nz(param.ep_face_node);
    
    % r4 eigenvector from \S3.5,
    % this is used to build flux for PII wave
    
    param.r4_p_ep_adj = zeros(num_pe_face_node, param.Nfields);
    
    param.r4_p_ep_adj(:, 1) = param.lambda_epM + 2*param.mu_fr_epM.*nx(param.ep_face_node).^2  + param.alpha_epM.*param.gamma2_epM.*param.M_epM;
    param.r4_p_ep_adj(:, 2) = param.lambda_epM + 2*param.mu_fr_epM.*ny(param.ep_face_node).^2  + param.alpha_epM.*param.gamma2_epM.*param.M_epM;
    param.r4_p_ep_adj(:, 3) = param.lambda_epM + 2*param.mu_fr_epM.*nz(param.ep_face_node).^2  + param.alpha_epM.*param.gamma2_epM.*param.M_epM;
    
    param.r4_p_ep_adj(:, 4) = 2*param.mu_fr_epM.*nx(param.ep_face_node).*ny(param.ep_face_node);
    param.r4_p_ep_adj(:, 5) = 2*param.mu_fr_epM.*ny(param.ep_face_node).*nz(param.ep_face_node);
    param.r4_p_ep_adj(:, 6) = 2*param.mu_fr_epM.*nx(param.ep_face_node).*nz(param.ep_face_node);
    
    param.r4_p_ep_adj(:, 7) = -param.alpha_epM.*param.M_epM - param.gamma2_epM.*param.M_epM;
    
    param.r4_p_ep_adj(:, 8) = param.cpII_i_epM.*nx(param.ep_face_node);
    param.r4_p_ep_adj(:, 9) = param.cpII_i_epM.*ny(param.ep_face_node);
    param.r4_p_ep_adj(:,10) = param.cpII_i_epM.*nz(param.ep_face_node);
    param.r4_p_ep_adj(:,11) = param.gamma2_epM.*param.cpII_i_epM.*nx(param.ep_face_node);
    param.r4_p_ep_adj(:,12) = param.gamma2_epM.*param.cpII_i_epM.*ny(param.ep_face_node);
    param.r4_p_ep_adj(:,13) = param.gamma2_epM.*param.cpII_i_epM.*nz(param.ep_face_node);

end % if num_interface_elt > 0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%direction cosines
nnx = param.pwx/sqrt(param.pwx^2+param.pwy^2+param.pwz^2); 
nny = param.pwy/sqrt(param.pwx^2+param.pwy^2+param.pwz^2); 
nnz = param.pwz/sqrt(param.pwx^2+param.pwy^2+param.pwz^2); % THIS OVERWRITES MATLAB nnz FUNCTION!

if num_elastic_elt > 0

    % used to build elastic plane wave in plane_wave_elastic
    Pi = nnx*A_e_adj(:,:,1) + nny*B_e_adj(:,:,1) + nnz*C_e_adj(:,:,1);
    [e_vec,e_val] = eig(Pi, Q_e_adj(:,:,1));
    e_val = diag(e_val);
    % sort eigenvalues and eigenvectors
    [~, ind] = sort(e_val);
    param.e_vec_plane_wave_ela_adj = e_vec(:,ind);
    param.c_plane_wave_ela_adj = e_val(ind);

    % Now inflate if nexessary by adding zero rows in positions 7, 10:13

    param.e_vec_plane_wave_e_adj = zeros(param.Nfields, param.Nfields_e);
    param.e_vec_plane_wave_e_adj(param.elastic_fields, :) = param.e_vec_plane_wave_ela_adj;
    
end % if num_elastic_elt > 0

% plane wave
if num_inviscid_elt > 0

    % used to build inviscid poroelatic plane wave in plane_wave_inviscid 
    Pi = nnx*param.A_p_adj(:,:,1) + nny*param.B_p_adj(:,:,1) + nnz*param.C_p_adj(:,:,1);
    [e_vec,e_val] = eig(Pi, param.Q_p(:,:,1));
    e_val = diag(e_val);
    % sort eigenvalues and eigenvectors
    [~, ind] = sort(e_val);
    param.e_vec_plane_wave_in_adj = e_vec(:,ind);
    param.c_plane_wave_in_adj = e_val(ind);

end % if num_inviscid_elt > 0
