function face_flux_ep = get_face_flux_ep(fieldM, fieldP, param)
% This code computes the downwind flux at the interface between elastic and
% poroelastic domains and returns the flux INTO the poroelastic domain.
% Note that get_face_flux_ep_pe.m returned both face_flux_ep and 
% face_flux_pe with the convention that the elastic elements were interior 
% point (minus) while the poroelastic elements were exterior (negative) 
% points. This presented a problem for parallelisation, so now both fluxes 
% are treated separately as fluxes INTO the respective interior elements. 
% This code treats the poroelastic element as an interior element and 
% elastic element as exterior and parallels the formal structure of
% get_face_flux_pp
% CONVENTIONS:
% param.ep_face_node contains all face node numbers where the
% interior node is on poroelastic element and the corresponding
% exterior node is on an elastic element
% face_flux_ep will hold flux from elastic to poroelastic elements
% 16 August 2018

persistent call_number;
if isempty(call_number)
    call_number = 1;
else
    call_number = call_number + 1;
end

% Return empty array if there are no pe faces

if isempty(param.pe_face_node)
    face_flux_ep = [];  
    return
end

Globals3D;

r1_p  = param.r1_p_ep;
r4_p  = param.r4_p_ep;

% Extract outward normals for poroelastic-elastic faces
onx = nx(param.ep_face_node);
ony = ny(param.ep_face_node);
onz = nz(param.ep_face_node);

% collect fields, poroelastic
e11M  = fieldM(param.ep_face_node,1);
e22M  = fieldM(param.ep_face_node,2);
e33M  = fieldM(param.ep_face_node,3);
e12M  = fieldM(param.ep_face_node,4);
e23M  = fieldM(param.ep_face_node,5);
e13M  = fieldM(param.ep_face_node,6);
zetaM = fieldM(param.ep_face_node,7);
vxM   = fieldM(param.ep_face_node,8);
vyM   = fieldM(param.ep_face_node,9);
vzM   = fieldM(param.ep_face_node,10);
vfxM  = fieldM(param.ep_face_node,11);
vfyM  = fieldM(param.ep_face_node,12);
vfzM  = fieldM(param.ep_face_node,13);

% P are elastic fields
e11P  = fieldP(param.ep_face_node,1);
e22P  = fieldP(param.ep_face_node,2);
e33P  = fieldP(param.ep_face_node,3);
e12P  = fieldP(param.ep_face_node,4);
e23P  = fieldP(param.ep_face_node,5);
e13P  = fieldP(param.ep_face_node,6);
% field 7 not used in elastic model
vxP   = fieldP(param.ep_face_node,8);
vyP   = fieldP(param.ep_face_node,9);
vzP   = fieldP(param.ep_face_node,10);

% compute velocity terms
% ndv = n*[vxM;vyM;vzM] - n*[vxP;vyP;vzP]
ndv  = ( onx.*vxM + ony.*vyM + onz.*vzM ) - ...
       ( onx.*vxP + ony.*vyP + onz.*vzP );
   
ndvf = onx.*vfxM + ony.*vfyM + onz.*vfzM;

% [[T,S]] = T^-n^- + S^+n^+
% ndTSn = n'[[T,S]]n
ndTSn = ...
( 2*param.mu_fr_epM.*( ...
   onx.*onx.*e11M + ...
   ony.*ony.*e22M + ...
   onz.*onz.*e33M + ...
 2*onx.*ony.*e12M + ...
 2*ony.*onz.*e23M + ...
 2*onx.*onz.*e13M   ...
) + ...
param.lambda_epM.*(e11M+e22M+e33M) - ...
param.alpha_epM.*param.M_epM.*zetaM ...       
) ...
- ... 
( 2*param.mu_epP.*( ...
  onx.*onx.*e11P + ...
  ony.*ony.*e22P + ...
  onz.*onz.*e33P + ...
2*onx.*ony.*e12P + ...
2*ony.*onz.*e23P + ...
2*onx.*onz.*e13P ...
) +...
param.lambda_epP.*(e11P+e22P+e33P) ...
);

% S-wave
TSnx =  (2*param.mu_fr_epM.*e11M + param.lambda_epM.*(e11M+e22M+e33M) - param.alpha_epM.*param.M_epM.*zetaM).*onx + ...
        2*param.mu_fr_epM.*e12M.*ony + ...
        2*param.mu_fr_epM.*e13M.*onz - ...
       (2*param.mu_epP.*e11P + param.lambda_epP.*(e11P+e22P+e33P)).*onx  -...
        2*param.mu_epP.*e12P.*ony - ...
        2*param.mu_epP.*e13P.*onz;  

TSny =   2*param.mu_fr_epM.*e12M.*onx + ...
        (2*param.mu_fr_epM.*e22M + param.lambda_epM.*(e11M+e22M+e33M) - param.alpha_epM.*param.M_epM.*zetaM).*ony + ...
         2*param.mu_fr_epM.*e23M.*onz - ...
         2*param.mu_epP.*e12P.*onx - ...
        (2*param.mu_epP.*e22P + param.lambda_epP.*(e11P+e22P+e33P)).*ony -...
         2*param.mu_epP.*e23P.*onz;

TSnz =  2*param.mu_fr_epM.*e13M.*onx +...
        2*param.mu_fr_epM.*e23M.*ony +...
       (2*param.mu_fr_epM.*e33M + param.lambda_epM.*(e11M+e22M+e33M) - param.alpha_epM.*param.M_epM.*zetaM).*onz -... 
        2*param.mu_epP.*e13P.*onx -...
        2*param.mu_epP.*e23P.*ony -...
       (2*param.mu_epP.*e33P + param.lambda_epP.*(e11P+e22P+e33P)).*onz;
        
      
% these are cross product terms nx(nx[[S]]) in \S3.5 
nnTSnx = onx.*(TSny.*ony + TSnz.*onz) - TSnx.*(ony.^2 + onz.^2); 

nnTSny = ony.*(TSnx.*onx + TSnz.*onz) - TSny.*(onx.^2 + onz.^2); 

nnTSnz = onz.*(TSnx.*onx + TSny.*ony) - TSnz.*(onx.^2 + ony.^2);

% this is nx(nx[v_s])
nndvx =  onx.*(ony.*vyM + onz.*vzM) - vxM.*(ony.^2 + onz.^2) -...
        (onx.*(ony.*vyP + onz.*vzP) - vxP.*(ony.^2 + onz.^2));
    
nndvy =  ony.*(onx.*vxM + onz.*vzM) - vyM.*(onx.^2 + onz.^2) -...
        (ony.*(onx.*vxP + onz.*vzP) - vyP.*(onx.^2 + onz.^2));
    
nndvz =  onz.*(onx.*vxM + ony.*vyM) - vzM.*(onx.^2 + ony.^2) -...
        (onz.*(onx.*vxP + ony.*vyP) - vzP.*(onx.^2 + ony.^2));


% storage for the elastic-to-poroelastic flux terms

face_flux_ep = zeros(numel(param.ep_face_node), param.Nfields);

tmp_i_S = -param.cs_i_epM./(param.cs_epP.*param.mu_fr_epM + param.cs_i_epM.*param.mu_epP);

face_flux_ep(:,1) = tmp_i_S.*onx.*(param.cs_epP.*nnTSnx + param.mu_epP.*nndvx);
face_flux_ep(:,2) = tmp_i_S.*ony.*(param.cs_epP.*nnTSny + param.mu_epP.*nndvy);
face_flux_ep(:,3) = tmp_i_S.*onz.*(param.cs_epP.*nnTSnz + param.mu_epP.*nndvz);
face_flux_ep(:,4) = 0.5*tmp_i_S.*(...
    param.cs_epP.*(ony.*nnTSnx + onx.*nnTSny) +...
    param.mu_epP.*(ony.*nndvx + onx.*nndvy) ...
    );
face_flux_ep(:,5) = 0.5*tmp_i_S.*(...
    param.cs_epP.*(onz.*nnTSny + ony.*nnTSnz) + ...
    param.mu_epP.*(onz.*nndvy + ony.*nndvz) ...
    );
face_flux_ep(:,6) = 0.5*tmp_i_S.*(...
    param.cs_epP.*(onz.*nnTSnx + onx.*nnTSnz) + ...
    param.mu_epP.*(onz.*nndvx + onx.*nndvz) ...
    );
%face_flux(:,7) = 0;
face_flux_ep(:,8)  =  tmp_i_S.*param.cs_i_epM.*( param.cs_epP.*nnTSnx + param.mu_epP.*nndvx );
face_flux_ep(:,9)  =  tmp_i_S.*param.cs_i_epM.*( param.cs_epP.*nnTSny + param.mu_epP.*nndvy );
face_flux_ep(:,10) =  tmp_i_S.*param.cs_i_epM.*( param.cs_epP.*nnTSnz + param.mu_epP.*nndvz );
face_flux_ep(:,11) = -tmp_i_S.*param.cs_i_epM.*param.rho_f_epM./param.m_epM.*( param.cs_epP.*nnTSnx + param.mu_epP.*nndvx );
face_flux_ep(:,12) = -tmp_i_S.*param.cs_i_epM.*param.rho_f_epM./param.m_epM.*( param.cs_epP.*nnTSny + param.mu_epP.*nndvy );
face_flux_ep(:,13) = -tmp_i_S.*param.cs_i_epM.*param.rho_f_epM./param.m_epM.*( param.cs_epP.*nnTSnz + param.mu_epP.*nndvz );

%fprintf('In get_face_flux_ep, step number %d\n', call_number);
%keyboard

% assemble fast P-wave flux
tmp_fastP = param.cpI_i_epM.*(param.d11_ep.*ndTSn + param.d12_ep.*ndv  + param.d13_ep.*ndvf);

for fld1 = 1:param.Nfields
    face_flux_ep(:,fld1) = face_flux_ep(:,fld1) + tmp_fastP.*r1_p(:, fld1);
end

% assemble slow P-wave flux
tmp_slowP = param.cpII_i_epM.*(param.d21_ep.*ndTSn + param.d22_ep.*ndv  + param.d23_ep.*ndvf);

for fld1 = 1:param.Nfields
    face_flux_ep(:,fld1) = face_flux_ep(:,fld1) + tmp_slowP.*r4_p(:, fld1);
end
%fprintf('Reached end of get_face_flux_ep, step number %d\n', call_number);
%keyboard