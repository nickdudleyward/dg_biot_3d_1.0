function param = setup_adjoint_wavefield_source(param)

Globals3D;

% Assemble adjoint wavefield source by taking the difference between the
% adjoint model solid velocities vx, vy, vz and the data velocity measurements

h5_file = fullfile(param.problem_data_dir, param.save_adjoint_filename);
tmp = h5read(h5_file, '/receiver_adj_model') - h5read(h5_file, '/receiver_data');
param.source_magnitude_adj_wavefield = tmp(:,:,8:10);

% Times at which receiver values were recorded during forward run of
% adjoint model in SETUP_ADJOINT_PROBLEM
% adjoint method run always starts at end of this run, so set
% param.initial_time accordingly
param.forward_step_time = h5read(h5_file, '/step_time');
param.initial_time = param.forward_step_time(end);

% For each receiver, assemble a vector akin to an inverse mass matrix: the
% Vandermonde matrix for the containing element times the transpose of the 
% Vandermonde matrix for the receiver point, transposed to store as a row of
% param.adjoint_source_IMM

param.adjoint_source_IMM = zeros(param.num_rec, Np);
for j=1:param.num_rec
    param.adjoint_source_IMM(j, :) = (1/J(1, param.elm_rec(j))*V*param.V_rec(j, :)')';
    % possibly cleaner version of previous line
    %param.adjoint_source_IMM(j, :) = param.V_rec(j, :)*V'/J(1, param.elm_rec(j));
end
