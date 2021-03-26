function face_flux_ee = get_face_flux_ee(fieldM, fieldP, param)
% Compute flux between adjacent elastic elements

persistent call_number;
if isempty(call_number)
    call_number = 1;
else
    call_number = call_number + 1;
end

% Return empty array if there are no ee faces

if isempty(param.ee_face_node)
    face_flux_ee = [];
    return
end

Globals3D;

% Extract outward normals for elastic-elastic faces

onx = nx(param.ee_face_node);
ony = ny(param.ee_face_node);
onz = nz(param.ee_face_node);

% TODO it might make more sense if we were passed 
% field[MP](param.ee_face_node,:) instead of the whole of field[MP]

% collect fields
e11M  = fieldM(param.ee_face_node, param.fld_e11);
e22M  = fieldM(param.ee_face_node, param.fld_e22);
e33M  = fieldM(param.ee_face_node, param.fld_e33);
e12M  = fieldM(param.ee_face_node, param.fld_e12);
e23M  = fieldM(param.ee_face_node, param.fld_e23);
e13M  = fieldM(param.ee_face_node, param.fld_e13);
vxM   = fieldM(param.ee_face_node, param.fld_vx);
vyM   = fieldM(param.ee_face_node, param.fld_vy);
vzM   = fieldM(param.ee_face_node, param.fld_vz);

e11P  = fieldP(param.ee_face_node, param.fld_e11);
e22P  = fieldP(param.ee_face_node, param.fld_e22);
e33P  = fieldP(param.ee_face_node, param.fld_e33);
e12P  = fieldP(param.ee_face_node, param.fld_e12);
e23P  = fieldP(param.ee_face_node, param.fld_e23);
e13P  = fieldP(param.ee_face_node, param.fld_e13);
vxP   = fieldP(param.ee_face_node, param.fld_vx);
vyP   = fieldP(param.ee_face_node, param.fld_vy);
vzP   = fieldP(param.ee_face_node, param.fld_vz);

% compute velocity terms
% ndv = n'[vxM;vyM;vzM] - n'[vxP;vyP;vzP]
ndv  = ( onx.*vxM + ony.*vyM + onz.*vzM ) - ...
       ( onx.*vxP + ony.*vyP + onz.*vzP );

% this is n'(S^- - S^+)n  
ndSn = ...
(   2*param.mu_eeM.*( ...
  onx.*onx.*e11M + ...
  ony.*ony.*e22M + ...
  onz.*onz.*e33M + ...
2*onx.*ony.*e12M + ...
2*ony.*onz.*e23M + ...
2*onx.*onz.*e13M ...
) +...
param.lambda_eeM.*(e11M+e22M+e33M) ...
) ...
- ...   
( 2*param.mu_eeP.*( ...
   onx.*onx.*e11P + ...
   ony.*ony.*e22P + ...
   onz.*onz.*e33P + ...
 2*onx.*ony.*e12P + ...
 2*ony.*onz.*e23P + ...
 2*onx.*onz.*e13P   ...
) + ...
param.lambda_eeP.*(e11P+e22P+e33P) ...       
);

% S-wave
Snx = (2*param.mu_eeM.*e11M + param.lambda_eeM.*(e11M+e22M+e33M)).*onx + ...
       2*param.mu_eeM.*e12M.*ony + ...
       2*param.mu_eeM.*e13M.*onz - ...
      (2*param.mu_eeP.*e11P + param.lambda_eeP.*(e11P+e22P+e33P)).*onx - ...
       2*param.mu_eeP.*e12P.*ony - ...
       2*param.mu_eeP.*e13P.*onz;

Sny = 2*param.mu_eeM.*e12M.*onx +...
     (2*param.mu_eeM.*e22M + param.lambda_eeM.*(e11M+e22M+e33M)).*ony +...
      2*param.mu_eeM.*e23M.*onz -...
      2*param.mu_eeP.*e12P.*onx -...
     (2*param.mu_eeP.*e22P + param.lambda_eeP.*(e11P+e22P+e33P)).*ony -...
      2*param.mu_eeP.*e23P.*onz;

Snz = 2*param.mu_eeM.*e13M.*onx +...
      2*param.mu_eeM.*e23M.*ony +...
     (2*param.mu_eeM.*e33M + param.lambda_eeM.*(e11M+e22M+e33M)).*onz -...
      2*param.mu_eeP.*e13P.*onx -...
      2*param.mu_eeP.*e23P.*ony -...
     (2*param.mu_eeP.*e33P + param.lambda_eeP.*(e11P+e22P+e33P)).*onz;

% these are cross product terms nx(nx[[S]]) in \S3.5 
nnSnx = onx.*(Sny.*ony + Snz.*onz) - Snx.*(ony.^2 + onz.^2); 

nnSny = ony.*(Snx.*onx + Snz.*onz) - Sny.*(onx.^2 + onz.^2); 

nnSnz = onz.*(Snx.*onx + Sny.*ony) - Snz.*(onx.^2 + ony.^2);

% this is nx(nx[v_s])
nndvx =  onx.*(ony.*vyM + onz.*vzM) - vxM.*(ony.^2 + onz.^2) -...
        (onx.*(ony.*vyP + onz.*vzP) - vxP.*(ony.^2 + onz.^2));
    
nndvy =  ony.*(onx.*vxM + onz.*vzM) - vyM.*(onx.^2 + onz.^2) -...
        (ony.*(onx.*vxP + onz.*vzP) - vyP.*(onx.^2 + onz.^2));
    
nndvz =  onz.*(onx.*vxM + ony.*vyM) - vzM.*(onx.^2 + ony.^2) -...
        (onz.*(onx.*vxP + ony.*vyP) - vzP.*(onx.^2 + ony.^2));

% storage for the flux terms.

face_flux_ee = zeros(numel(param.ee_face_node), param.Nfields);

% assemble S-wave flux
tmpS = -param.cs_eeM./(param.cs_eeP.*param.mu_eeM + param.cs_eeM.*param.mu_eeP);

face_flux_ee(:, param.fld_e11) = tmpS.*onx.*(param.cs_eeP.*nnSnx + param.mu_eeP.*nndvx);
face_flux_ee(:, param.fld_e22) = tmpS.*ony.*(param.cs_eeP.*nnSny + param.mu_eeP.*nndvy);
face_flux_ee(:, param.fld_e33) = tmpS.*onz.*(param.cs_eeP.*nnSnz + param.mu_eeP.*nndvz);
face_flux_ee(:, param.fld_e12) = 0.5*tmpS.*(...
    param.cs_eeP.*(ony.*nnSnx + onx.*nnSny) +...
    param.mu_eeP.*(ony.*nndvx + onx.*nndvy) ...
    );
face_flux_ee(:, param.fld_e23) = 0.5*tmpS.*(...
    param.cs_eeP.*(onz.*nnSny + ony.*nnSnz) +...
    param.mu_eeP.*(onz.*nndvy + ony.*nndvz) ...
    );
face_flux_ee(:, param.fld_e13) = 0.5*tmpS.*(...
    param.cs_eeP.*(onz.*nnSnx + onx.*nnSnz) +...
    param.mu_eeP.*(onz.*nndvx + onx.*nndvz) ...
    );
face_flux_ee(:, param.fld_vx) = tmpS.*param.cs_eeM.*(param.cs_eeP.*nnSnx + param.mu_eeP.*nndvx);
face_flux_ee(:, param.fld_vy) = tmpS.*param.cs_eeM.*(param.cs_eeP.*nnSny + param.mu_eeP.*nndvy);
face_flux_ee(:, param.fld_vz) = tmpS.*param.cs_eeM.*(param.cs_eeP.*nnSnz + param.mu_eeP.*nndvz);

% assemble P-wave flux
tmpP = param.cp_eeM.*(param.d11_ee.*ndSn + param.d12_ee.*ndv);

for fld1 = param.elastic_fields
    face_flux_ee(:,fld1) = face_flux_ee(:,fld1) + tmpP.*param.r1_e_ee(:, fld1);
end

% fprintf('Reached end of face_flux_ee, step number %d\n', call_number);
% keyboard