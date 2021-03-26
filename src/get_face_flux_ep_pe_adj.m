function [face_flux_ep, face_flux_pe] = get_face_flux_ep_pe_adj(fieldM, fieldP, param)
% Compute fluxes in both directions between adjacent elastic and 
% poroelastic elements. Here elastic elements are interior and poroelastic
% elements are exterior. The assembly of the fluxes follows the paper:
% for the elastic downwind flux:
%
% (Pi q)^* = (Pi^-)q^- + (sum weight_i * r_i^-) 
%
% where the eigenvectors correspond to incoming (-) elastic eigenvectors
% For the poroelastic downwind flux
% 
% (Pi q)^* = (Pi^+)q^+ - (sum weight_i * r_i^+) 
%
% where the eigenvectors correspond to incoming (+) poroelatic eigenvectors
% The weights param.dij_pe are constructed using the conventions in the
% paper, \S5.1
% This code caused a great deal of inconvenience, the main culprit
% was ndvf = onx.*vfxP + ony.*vfyP + onz.*vfzP;
% before we left in -: ndvf = -(onx.*vfxP + ony.*vfyP + onz.*vfzP);
% this gets removed in definition of flux, see paper
% CONVENTIONS:
% param.pe_face_node contains all face node numbers where the
% interior node is on an elastic element and the corresponding
% exterior node is on a poroelastic element
% face_flux_ep will hold flux from elastic to poroelastic elements
% face_flux_pe will hold flux from poroelastic to elastic elements
% 16 August 2018

% 6 feb 2019
% implemented for adjoint formulation

persistent call_number;
if isempty(call_number)
    call_number = 1;
else
    call_number = call_number + 1;
end

% Return empty array if there are no pe faces

if isempty(param.pe_face_node)
    face_flux_ep = [];
    face_flux_pe = [];
    return
end

Globals3D;

r1_ie_adj = param.r1_ie_pe_adj;
r1_p_adj  = param.r1_p_pe_adj;
r4_p_adj  = param.r4_p_pe_adj;
% WHY?????? presumable because of normal conventions
% r1_p_adj  = param.r1_p_ep_adj;
% r4_p_adj  = param.r4_p_ep_adj;

% Extract matrices of outward normals from elastic to poroelastic 
% neighbour. 

onx = nx(param.pe_face_node);
ony = ny(param.pe_face_node);
onz = nz(param.pe_face_node);

% Work initially with poro-elastic-to-elasic fluxes. onx, ony, onz
% get changed in sign when we assemble the elastic-to-poroelastic fluxes

% In this section the minus M corresponds to an elastic field and P to the
% adjacent poroelastic field
% TODO it might make more sense if we were passed 
% field[MP](param.pe_face_node,:) instead of field[MP] in its
% entirety

% collect fields, elastic
s11M  = fieldM(param.pe_face_node,1);
s22M  = fieldM(param.pe_face_node,2);
s33M  = fieldM(param.pe_face_node,3);
s12M  = fieldM(param.pe_face_node,4);
s23M  = fieldM(param.pe_face_node,5);
s13M  = fieldM(param.pe_face_node,6);
% field 7 not used in elastic model
vxM   = fieldM(param.pe_face_node,8);
vyM   = fieldM(param.pe_face_node,9);
vzM   = fieldM(param.pe_face_node,10);

% P are poroelastic fields
s11P  = fieldP(param.pe_face_node,1);
s22P  = fieldP(param.pe_face_node,2);
s33P  = fieldP(param.pe_face_node,3);
s12P  = fieldP(param.pe_face_node,4);
s23P  = fieldP(param.pe_face_node,5);
s13P  = fieldP(param.pe_face_node,6);
%pP    = fieldP(param.pe_face_node,7); NOT USED
vxP   = fieldP(param.pe_face_node,8);
vyP   = fieldP(param.pe_face_node,9);
vzP   = fieldP(param.pe_face_node,10);
vfxP  = fieldP(param.pe_face_node,11);
vfyP  = fieldP(param.pe_face_node,12);
vfzP  = fieldP(param.pe_face_node,13);

% compute velocity terms
% ndv = n*[vxM;vyM;vzM] - n*[vxP;vyP;vzP]
ndv  = ( onx.*vxM + ony.*vyM + onz.*vzM ) - ...
       ( onx.*vxP + ony.*vyP + onz.*vzP );

ndvf = onx.*vfxP + ony.*vfyP + onz.*vfzP;

% [[S,T]] = S^-n^- + T^+n^+
% ndSTn = n'[[S,T]]n
% ndSTn = ...
% (   2*param.mu_peM.*( ...
%   onx.*onx.*e11M + ...
%   ony.*ony.*e22M + ...
%   onz.*onz.*e33M + ...
% 2*onx.*ony.*e12M + ...
% 2*ony.*onz.*e23M + ...
% 2*onx.*onz.*e13M ...
% ) +...
% param.lambda_peM.*(e11M+e22M+e33M) ...
% ) ...
% - ...   
% ( 2*param.mu_fr_peP.*( ...
%    onx.*onx.*e11P + ...
%    ony.*ony.*e22P + ...
%    onz.*onz.*e33P + ...
%  2*onx.*ony.*e12P + ...
%  2*ony.*onz.*e23P + ...
%  2*onx.*onz.*e13P   ...
% ) + ...
% param.lambda_peP.*(e11P+e22P+e33P) - ...
% param.alpha_peP.*param.M_peP.*zetaP ...       
% );

ndSTn = ...
( ...
  onx.*onx.*s11M + ...
  ony.*ony.*s22M + ...
  onz.*onz.*s33M + ...
2*onx.*ony.*s12M + ...
2*ony.*onz.*s23M + ...
2*onx.*onz.*s13M  ...
) ...
- ...   
( ...
   onx.*onx.*s11P + ...
   ony.*ony.*s22P + ...
   onz.*onz.*s33P + ...
 2*onx.*ony.*s12P + ...
 2*ony.*onz.*s23P + ...
 2*onx.*onz.*s13P ...     
);

% S-wave
% STnx = (2*param.mu_peM.*e11M + param.lambda_peM.*(e11M+e22M+e33M)).*onx + ...
%         2*param.mu_peM.*e12M.*ony + ...
%         2*param.mu_peM.*e13M.*onz - ...  
%        (2*param.mu_fr_peP.*e11P + param.lambda_peP.*(e11P+e22P+e33P) - param.alpha_peP.*param.M_peP.*zetaP).*onx - ...
%         2*param.mu_fr_peP.*e12P.*ony - ...
%         2*param.mu_fr_peP.*e13P.*onz;  
%    
% STny =  2*param.mu_peM.*e12M.*onx +...
%        (2*param.mu_peM.*e22M + param.lambda_peM.*(e11M+e22M+e33M)).*ony +...
%         2*param.mu_peM.*e23M.*onz -...
%         2*param.mu_fr_peP.*e12P.*onx -...
%        (2*param.mu_fr_peP.*e22P + param.lambda_peP.*(e11P+e22P+e33P) - param.alpha_peP.*param.M_peP.*zetaP).*ony -...
%         2*param.mu_fr_peP.*e23P.*onz;
% 
% STnz =  2*param.mu_peM.*e13M.*onx +...
%         2*param.mu_peM.*e23M.*ony +...
%        (2*param.mu_peM.*e33M + param.lambda_peM.*(e11M+e22M+e33M)).*onz -...
%         2*param.mu_fr_peP.*e13P.*onx -...
%         2*param.mu_fr_peP.*e23P.*ony -...
%        (2*param.mu_fr_peP.*e33P + param.lambda_peP.*(e11P+e22P+e33P) - param.alpha_peP.*param.M_peP.*zetaP).*onz; 

STnx =  s11M.*onx + s12M.*ony + s13M.*onz -...
       (s11P.*onx + s12P.*ony + s13P.*onz);


STny =  s12M.*onx + s22M.*ony + s23M.*onz -...
       (s12P.*onx + s22P.*ony + s23P.*onz);

STnz =  s13M.*onx + s23M.*ony + s33M.*onz -...
       (s13P.*onx + s23P.*ony + s33P.*onz);
      
% these are cross product terms nx(nx[[S]]) in \S3.5 
nnSTnx = onx.*(STny.*ony + STnz.*onz) - STnx.*(ony.^2 + onz.^2); 

nnSTny = ony.*(STnx.*onx + STnz.*onz) - STny.*(onx.^2 + onz.^2); 

nnSTnz = onz.*(STnx.*onx + STny.*ony) - STnz.*(onx.^2 + ony.^2);

% this is nx(nx[v_s])
nndvx =  onx.*(ony.*vyM + onz.*vzM) - vxM.*(ony.^2 + onz.^2) -...
        (onx.*(ony.*vyP + onz.*vzP) - vxP.*(ony.^2 + onz.^2));
    
nndvy =  ony.*(onx.*vxM + onz.*vzM) - vyM.*(onx.^2 + onz.^2) -...
        (ony.*(onx.*vxP + onz.*vzP) - vyP.*(onx.^2 + onz.^2));
    
nndvz =  onz.*(onx.*vxM + ony.*vyM) - vzM.*(onx.^2 + ony.^2) -...
        (onz.*(onx.*vxP + ony.*vyP) - vzP.*(onx.^2 + ony.^2));

% Assembly of fluxes, \S5.2. Fisrt poroelastic-to-elastic flux terms then
% elastic-to-poroelastic flux terms. The first uses eignvectors for the
% elastic element while the second for the poroelastic element

% BEGIN ASSEMBLY POROELASTIC-TO-ELASTIC WAVE FLUX
% This parallels assembly in get_face_flux_ee_adj.m
% storage for poroelastic-to-elastic flux terms


%--------------------------------------------------------------------------

face_flux_pe = zeros(numel(param.pe_face_node), param.Nfields);

% Elastic wave
% assemble S-wave flux
tmpS = -param.cs_peM./(param.cs_i_peP.*param.mu_peM + param.cs_peM.*param.mu_fr_peP);

% think i missed this out
tmp2 = 2*param.mu_peM;

face_flux_pe(:,1) = tmp2.*tmpS.*onx.*(param.cs_i_peP.*nnSTnx + param.mu_fr_peP.*nndvx);
face_flux_pe(:,2) = tmp2.*tmpS.*ony.*(param.cs_i_peP.*nnSTny + param.mu_fr_peP.*nndvy);
face_flux_pe(:,3) = tmp2.*tmpS.*onz.*(param.cs_i_peP.*nnSTnz + param.mu_fr_peP.*nndvz);
face_flux_pe(:,4) = 0.5*tmp2.*tmpS.*(...
     param.cs_i_peP.*(ony.*nnSTnx + onx.*nnSTny) +...
     param.mu_fr_peP.*(ony.*nndvx + onx.*nndvy) ...
    );
face_flux_pe(:,5) = 0.5*tmp2.*tmpS.*(...
     param.cs_i_peP.*(onz.*nnSTny + ony.*nnSTnz) +...
     param.mu_fr_peP.*(onz.*nndvy + ony.*nndvz) ...
    );
face_flux_pe(:,6) = 0.5*tmp2.*tmpS.*(...
     param.cs_i_peP.*(onz.*nnSTnx + onx.*nnSTnz) +...
     param.mu_fr_peP.*(onz.*nndvx + onx.*nndvz) ...
    );
% row 7 not used, stays as initialised to zero
face_flux_pe(:,8)  = tmpS.*param.cs_peM.*(param.cs_i_peP.*nnSTnx + param.mu_fr_peP.*nndvx);
face_flux_pe(:,9)  = tmpS.*param.cs_peM.*(param.cs_i_peP.*nnSTny + param.mu_fr_peP.*nndvy);
face_flux_pe(:,10) = tmpS.*param.cs_peM.*(param.cs_i_peP.*nnSTnz + param.mu_fr_peP.*nndvz);

% assemble elastic P-wave flux

tmpP = param.cp_peM.*(param.d11_pe.*ndSTn + param.d12_pe.*ndv + param.d13_pe.*ndvf);

for fld1 = 1:param.Nfields
    face_flux_pe(:,fld1) = face_flux_pe(:,fld1) + tmpP.*r1_ie_adj(:, fld1);
end
% END ASSEMBLY POROELASTIC-TO-ELASTIC WAVE FLUX

% BEGIN ASSEMBLY ELASTIC-TO-POROELASTIC WAVE FLUX
% This needs assembling with great care since working with exterior element

% if the pe and ep face node lists are not sorted in the correct order,
% there's no point trying this!

if ~param.align_ep_pe_nodes
    face_flux_ep = [];
    return
end

% storage for the elastic-to-poroelastic flux terms

face_flux_ep = zeros(numel(param.pe_face_node), param.Nfields);

tmp_i_S = -param.cs_i_peP./(param.cs_i_peP.*param.mu_peM + param.cs_peM.*param.mu_fr_peP);

% think i missed this out
% NEED TO CHECK THIS MIGHT NEED tmp2 = 2*param.mu_fr something P CHECK
% DERIVATION
tmp2 = 2*param.mu_fr_peP;

% S-wave flux
% this flux is constructed using the conventions of the paper, in
% particular using inward pointing normals to poroelastic elements

face_flux_ep(:,1) = tmp2.*tmp_i_S.*onx.*(param.cs_peM.*nnSTnx - param.mu_peM.*nndvx);
face_flux_ep(:,2) = tmp2.*tmp_i_S.*ony.*(param.cs_peM.*nnSTny - param.mu_peM.*nndvy);
face_flux_ep(:,3) = tmp2.*tmp_i_S.*onz.*(param.cs_peM.*nnSTnz - param.mu_peM.*nndvz);
face_flux_ep(:,4) = 0.5*tmp2.*tmp_i_S.*(...
    param.cs_peM.*(ony.*nnSTnx + onx.*nnSTny) -...
    param.mu_peM.*(ony.*nndvx + onx.*nndvy) ...
    );
face_flux_ep(:,5) = 0.5*tmp2.*tmp_i_S.*(...
    param.cs_peM.*(onz.*nnSTny + ony.*nnSTnz) - ...
    param.mu_peM.*(onz.*nndvy + ony.*nndvz) ...
    );
face_flux_ep(:,6) = 0.5*tmp2.*tmp_i_S.*(...
    param.cs_peM.*(onz.*nnSTnx + onx.*nnSTnz) - ...
    param.mu_peM.*(onz.*nndvx + onx.*nndvz) ...
    );
%face_flux(:,7) = 0;
face_flux_ep(:,8)  = -tmp_i_S.*param.cs_i_peP.*( param.cs_peM.*nnSTnx - param.mu_peM.*nndvx );
face_flux_ep(:,9)  = -tmp_i_S.*param.cs_i_peP.*( param.cs_peM.*nnSTny - param.mu_peM.*nndvy );
face_flux_ep(:,10) = -tmp_i_S.*param.cs_i_peP.*( param.cs_peM.*nnSTnz - param.mu_peM.*nndvz );
face_flux_ep(:,11) =  tmp_i_S.*param.cs_i_peP.*param.rho_f_peP./param.m_peP.*( param.cs_peM.*nnSTnx - param.mu_peM.*nndvx );
face_flux_ep(:,12) =  tmp_i_S.*param.cs_i_peP.*param.rho_f_peP./param.m_peP.*( param.cs_peM.*nnSTny - param.mu_peM.*nndvy );
face_flux_ep(:,13) =  tmp_i_S.*param.cs_i_peP.*param.rho_f_peP./param.m_peP.*( param.cs_peM.*nnSTnz - param.mu_peM.*nndvz );

% assemble fast P-wave flux
tmp_fastP = param.cpI_i_peP.*(param.d31_pe.*ndSTn + param.d32_pe.*ndv  + param.d33_pe.*ndvf);

for fld1 = 1:param.Nfields
    face_flux_ep(:,fld1) = face_flux_ep(:,fld1) + tmp_fastP.*r1_p_adj(:, fld1);
end

% assemble slow P-wave flux
tmp_slowP = param.cpII_i_peP.*(param.d21_pe.*ndSTn + param.d22_pe.*ndv  + param.d23_pe.*ndvf);

for fld1 = 1:param.Nfields
    face_flux_ep(:,fld1) = face_flux_ep(:,fld1) + tmp_slowP.*r4_p_adj(:, fld1);
end

% see (154) paper for sign change CHECK LATER

face_flux_ep = - face_flux_ep;

%fprintf('Reached end of get_face_flux_ep_pe, step number %d\n', call_number);
%keyboard
% END ASSEMBLY ELASTIC-TO-POROELASTIC WAVE FLUX