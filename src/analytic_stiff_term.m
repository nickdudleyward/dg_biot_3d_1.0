function field = analytic_stiff_term(field, param, dt)
fprintf('Analytic stiff term\n');
% This is all terribly inefficient but it's written in the hope of being
% readable and the costs of running this are insignificant compared to the
% main RHS code.

% If there are no LF elements, immediately return empty

if isempty(field)
    return
end

%  extract the last six fields, i.e. the velocities
vx  = field(:, :, 8);
vy  = field(:, :, 9);
vz  = field(:, :, 10);
vfx = field(:, :, 11);
vfy = field(:, :, 12);
vfz = field(:, :, 13);

% These are the two (apart from 0 and 1) terms in the exponential of the 
% matrix representing the stiff part of the equation

t1 = exp(param.stiff_eigenval*dt);
t2 = -(t1-1).*param.rho_f_lf./param.rho_a_lf;

% This is the matrix that we want to multiply on the the velocity fields
% See Maple worksheet stiff_terms.mw

% [[1, 0, 0, t2, 0, 0];
%  [0, 1, 0, 0, t2, 0]; 
%  [0, 0, 1, 0, 0, t2]; 
%  [0, 0, 0, t1, 0, 0]; 
%  [0, 0, 0, 0, t1, 0];
%  [0, 0, 0, 0, 0, t1]]

%keyboard

% Note on transpose and .*
% t1 and t2 are column vectors, so t1' and t2' as used below are row
% vectors. When we .* a row vector onto a matrix, it multiplies the vector
% pointwise onto each row. Thus, t1 and t2 are calculated at element level
% and applied to each node on each element

vx_out = vx + t2'.*vfx;
vy_out = vy + t2'.*vfy;
vz_out = vz + t2'.*vfz;

vfx_out = t1'.*vfx;
vfy_out = t1'.*vfy;
vfz_out = t1'.*vfz;

field(:, :, 8)  = vx_out;
field(:, :, 9)  = vy_out;
field(:, :, 10) = vz_out;
field(:, :, 11) = vfx_out;
field(:, :, 12) = vfy_out;
field(:, :, 13) = vfz_out;
