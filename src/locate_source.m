function param = locate_source(param)

Globals3D;

% find element containing source and barycentric coordinates
% of source within element
[param.elm_sou, bary_sou] = tsearchn([VX',VY',VZ'], EToV, param.g_sou);
% Convert barycentrics to r, s, t coordinate system: H&W p.409
rst_sou = [-1,1,-1,-1; -1,-1,1,-1;-1,-1,-1,1]*bary_sou';
% Leaving this for the time being but if we move to multiple sources the
% following line will give columns containing r, s, t coords instead of 
% rows as in line above.
% This distinction caused problems in setup_receiver.
% SPE 2018/12/11
%rst_sou = bary_rec*[[-1, -1, -1]; [1, -1, -1]; [-1, 1, -1]; [-1, -1, 1]];
param.r_sou = rst_sou(1);
param.s_sou = rst_sou(2);
param.t_sou = rst_sou(3);

% Find model type (elastic, HF, LF, ...) of containing element

param.source_elt_model = param.elt_to_model_type(param.elm_sou);

% Find location of element within the list of elastic and within the
% list of poroelastic elements (needed because mass matrices are held 
% separately in param.Qi_ie and param.Qi_p)

% One of these will return a single index number, the other an empty vector

param.por_elt_sou = find(param.poroelastic_elt == param.elm_sou);
param.ela_elt_sou = find(param.elastic_elt == param.elm_sou);

