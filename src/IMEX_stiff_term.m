function res_field_IM = IMEX_stiff_term(r, res_field_EX, param)

%fprintf('IMEX stiff term\n');
% This function returns as res_field_IM
% (I-rA)^(-1)*A*y
% where A is the potentially stiff linear low-frequency dissipation term
% in the Biot equation, r is a scalar parameter (determined by the caller
% from the IMEX coefficients and time step) and y is received as parameter
% res_field_EX
%
% Only fields on the low-frequency element nodes are passed and returned

% Initialise; those field locations not affected by LF dissipation so 
% just hold zero, those affected are overwritten later

res_field_IM = zeros(size(res_field_EX));

% If there are no LF elements, return empty matrix at this point

if isempty(res_field_IM)
    return
end

% The matrix we want to multiply by is zero apart from this submatrix in 
% rows 8:13 (vx, vy, vz, vfx, vfy, vfz) and columns 11:13 (vfx, vfy, vfz)
% (see Maple worksheet stiff_terms.mw)
%
% eta/(r*eta*rho_a+k*m*rho_a-k*rho_f^2) *
% [rho_f,   0,      0      ]
% [0,       rho_f,  0      ]
% [0,       0,      rho_f  ]
% [-rho_a,  0,      0      ]
% [0,       -rho_a, 0      ]
% [0,       0,      -rho_a ]

% so there are actually only two non-zero terms present; call these ts
% in the upper half (affecting solid velocities) and tf in the lower
% half (affecting fluid velocities). Calculate these for each element 
% (they can't be cached because the parameter r varies).
% Transposing ts and tf makes them row vectors, which is important below.

tmp = param.eta_lf./( ...
    r*param.eta_lf.*param.rho_a_lf + ...
    param.k_lf.*param.m_lf.*param.rho_a_lf - ...
    param.k_lf.*param.rho_f_lf.^2 ...
);
ts = ( tmp.*param.rho_f_lf)'; % upper block, multiply onto solid velocities
tf = (-tmp.*param.rho_a_lf)'; % lower block, multiply onto fluid velocities

% Now apply multiplications. Note on .*
% ts and tf are row vectors. When we .* a row vector onto a matrix, it 
% multiplies the vector pointwise onto each row. Thus, ts and tf are 
% calculated at element level and applied to each node on each element

res_field_IM(:,:, 8) = ts.*res_field_EX(:,:,11); % vx
res_field_IM(:,:, 9) = ts.*res_field_EX(:,:,12); % vy
res_field_IM(:,:,10) = ts.*res_field_EX(:,:,13); % vz

res_field_IM(:,:,11) = tf.*res_field_EX(:,:,11); % vfx
res_field_IM(:,:,12) = tf.*res_field_EX(:,:,12); % vfy
res_field_IM(:,:,13) = tf.*res_field_EX(:,:,13); % vfz
