function param = get_inviscid_wavespeeds(param)
% computes inviscid wave speeds based on material parameters defined in
% param. these are used in the numerical schemes for all cases (inviscid, low/high freq dissipative)

% flux stuff and wave speeds.  Wave speeds for the non-dissipative case
% follow conventions in paper
Z1 = (param.m.*param.rho_a - param.rho_f.^2).^(-1);
Z2 = -2.*param.rho_f.*param.alpha.*param.M + param.rho_a.*param.M + param.m.*param.lambda + 2.*param.m.*param.mu_fr;
Z3 = param.rho_a.*(4.*param.alpha.^2.*param.m-4.*param.alpha.*param.rho_f+param.rho_a).*param.M.^2-...
    2.*(2.*param.alpha.*param.m.*param.rho_f+param.m.*param.rho_a-2.*param.rho_f.^2).*param.M.*(2.*param.mu_fr+param.lambda)+...
    param.m.^2.*(2.*param.mu_fr+param.lambda).^2;
Z4 = param.rho_a.*param.M - param.m.*param.lambda-2.*param.m.*param.mu_fr;
Z5 = 2.*(param.alpha.*param.m-param.rho_f).*param.M;
% inviscid/non-dissipative wave speeds
param.cpI_i  = sqrt(Z1/2.*(Z2+sqrt(Z3)));
param.cpII_i = sqrt(Z1/2.*(Z2-sqrt(Z3)));
param.cs_i   = sqrt(param.m.*Z1.*param.mu_fr);
% THESE ARE IMPORTANT
% gammas
param.gamma1 = (Z4+sqrt(Z3))./Z5;
param.gamma2 = (Z4-sqrt(Z3))./Z5;
