function driver_biot_3d(varargin)

% Problem file is normally passed as an argument. If no filename is passed 
% as an argument, use default below. Intended to be used interactively in
% the IDE
%default_problem_file = 'problem_adjoint_forward_solver_test';
default_problem_file = 'problem_adjoint_method_test_2';
%default_problem_file = 'problem_HFCT';

% This overrides any filename passed as an argument and cuts out any 
% checking of number of arguments. Intended to be used with great care,
% if at all.

override_problem_file = '';

% Find name of problem-specific config file

if ~isempty(override_problem_file)
    problem_file = override_problem_file;
elseif nargin == 0
    if ~isempty(default_problem_file)
        problem_file = default_problem_file;
    else
        error('No problem file specified on command line or in default_problem_file')
    end
elseif nargin == 1
    problem_file = varargin{1};
else
    error('More than one argument specified')
end

% and turn filename into a function handle

problem = str2func(problem_file);

% find path of parent directory of this script's location (each fileparts 
% removes one name from the end of the directory of the executing script).

parent_path = fileparts(fileparts(mfilename('fullpath')));

% and append parent/lib/hest_warb, parent/src and parent/examples
% to the path

addpath(fullfile(parent_path, 'lib/hest_warb'));
addpath(fullfile(parent_path, 'src'));
addpath(fullfile(parent_path, 'examples'));

%keyboard

% Load integer constants for magic numbers (to make closer parallel
% between Matlab and C code). This effectively initialise the param
% struct.

param = setup_magic();

% Begin first-stage user configuration. This takes place before H&W's
% Setup3D.

% Set defaults for some param fields, before problem file is run
% Changes here will potentially affect all problem files!
% These are all overridden in the problem file where required
% (taking care to use the correct configuration stage)
                                  
% number of elastic fields and number of HF memory variables. Total
% number of fields set later (after domain setup) to 9 or 13 for
% elastic or poroelastic model

param.num_mem_hf = 3;
param.Nfields_e  = 9;

% Time-stepping configuration

param.initial_time = 0;           % Initial time
param.freq_sou = 2000;            % Source frequency
param.time_width = 1/param.freq_sou; % One wave cycle
param.num_steps = 0;              % If >0, override param.time_width and 
                                  % use param.num_steps time steps
param.cfl = 0.4;                  % timestep multiplier in CFL condition
param.steps_per_report = 1;       % print out message after each RK step

% Operator splitting and Runge-Kutta method (including IMEX;
% param.RK_method is a group name in RK.h5 containing weights etc.)
                                  
param.splitting = param.SPL_NONE; % no operator splitting
param.RK_method = 'CK4';          % Carpenter-Kennedy order 4

% Booleans controlling the type of simulation run

% If false, runs setup (to save for C framework) but no simulation
param.run_simulation = true;
% If above is true:
param.run_forward = true;         % forward problem solver, velocity/strain
param.run_forward_adj = false;    % forward problem solver, velocity/stress
param.run_adjoint_method = false; % adjoint method implementation

% Directories for loading and saving data
                                  
param.problem_data_dir = pwd;     % problem data files (HDF5)
param.mesh_dir = pwd;             % mesh data

% Filenames for loading and saving data

% for save and restart
param.save_state_filename = 'biot_3d_state.h5';
% step times, field and receiver data for adjoint method
param.save_adjoint_filename = 'biot_3d_adjoint_rec_data.h5';
% problem-specific data for parallel C framework
param.C_extra_filename = 'biot_3d_extra.h5';
% boundary type data for parallel C framework
param.bctype_filename =  'biot_3d_bctype.txt';
% saved field and receiver data for comparison with C framework output
param.output_data_filename = 'biot_3d_m.h5';

% Switches controlling (mostly) saving data in above files

param.save_for_restart = false;   % save final state for later restart
param.restart_from_saved = false; % load initial state from earlier save
% Control saving of final state data for later adjoint method run: can be 
% param.SAVE_RECEIVER_DATA or param.SAVE_ADJ_MODEL
param.save_adjoint_data = param.SAVE_ADJOINT_NONE;
% Control saving of final state data for C framework to read
param.save_params = false;
param.save_output_field = true;
param.save_output_rec = true;

% Initial states; empty is often better than undefined

% no receivers
param.g_rec = [];
% Correct setting for no HF elements; if HF elements exist, these must be
% initialised in problem file
param.mem_hf = [];
param.mem_hf_adj = [];
% field and adjoint field empty, so solver can be run without both fields
% being otherwise initialised
param.field = [];
param.field_adj = [];
% pure_elastic empty means infer from (non)existence of poro elements
param.pure_elastic = [];

% Other config parameters

param.align_ep_pe_nodes = false;  % If true, switches on earlier code
                                  % version that can't be paralleled
param.CT_model_type = param.NONE; % Model type, for convergence tests in
                                  % C framework only
param.weight =  [1 1 1 1];        % For reweighting parts of the plane wave
                                  % analytic solutions

% Run first-stage problem-specific configuration

param = problem(param.PRE_MESH_LOAD, param);

% Throw an error if any of these fields hasn't been assigned in param

check_mandatory(param, {'mesh', 'order'});

% set up angular frequency, from source frequency

param.omega_sou = 2*pi*param.freq_sou;

% default to having no sensitivity kernel calculations

param.num_kernels = 0;
param.sens_ker = [];

% End of first-stage problem configuration. Now set up H&W framework.

% Hesthaven-Warburton global variables
% clear global is important here for interactive use, otherwise values
% persist between runs: in particular, BCType holds its values from
% the previous session, causing an error.

clear global;
Globals3D;
N = param.order; % H&W global name

% Read mesh data. If mesh in the problem file specifies a .mat file, 
% this is assumed to contain (at least) EToV, VX, VY, VZ

% If param.mesh is the name of a .mat file, assume it has variables
% EToV, VX, VY, VZ and load these into the global variables of the same
% names. Otherwise, assume it's the basename of plain text files 
% *_element.txt and *_vertex.txt, containing the same data.

% Similarly, load domain data if present; this isn't needed yet,
% but it seems sensible to load all mesh data at the same time.

% TODO: if the mesh generator is also to supply boundary condition
% data (currently this is all done in problem config stage 3), it should be 
% loaded here. Note that, as currently coded, this will trip an error:
% the system requires BCType to be zero or empty before problem config 
% stage 3, when all the mapX and vmapX vectors should be set; it then
% reverse-engineers BCType and (if param.save_params is true) saves it 
% for the C framework to find.

if endsWith(param.mesh, '.mat')
    mesh_data = load(fullfile(param.mesh_dir, param.mesh));
    EToV = mesh_data{'EToV'};
    VX =  mesh_data{'VX'};
    VY =  mesh_data{'VY'};
    VZ =  mesh_data{'VZ'};
    if isfield(mesh_data, 'Domainindex')
        param.elt_to_domain = mesh_data{'Domainindex'};
    end
else
    EToV = dlmread(fullfile(param.mesh_dir, strcat(param.mesh, '_element.txt', ' ')));
    VV =   dlmread(fullfile(param.mesh_dir, strcat(param.mesh, '_vertex.txt',  ' ')));
    VX =   VV(:,1)';
    VY =   VV(:,2)';
    VZ =   VV(:,3)';
    try
        param.elt_to_domain = dlmread(fullfile(param.mesh_dir, strcat(param.mesh, '_domain.txt', ' ')));
    catch
    end
end

% Set H&W global variables to reflect mesh size ...

K = size(EToV, 1); % number of elements
Nv = size(VX, 2);  % number of vertices (not nodes)

% ... and set up Hesthaven-Warburton framework (which must happen after 
% mesh has been loaded)

StartUp3D;

%keyboard

% Continue setup prior to problem config stage 2

% Find minimal edge length across mesh, used for CFL condition. 
% TODO This should probably be a separate function; should h_min be saved 
% in param?; is this the right place to find this?

h2_min = 0;
for t=1:K
    E = EToV(t,:);
    % vector - transpose gives a matrix of pairwise differences, so h2
    % is 4x4 symmetric with zero diagonal
    h2 = (VX(E)-VX(E)').^2+(VY(E)-VY(E)').^2+(VZ(E)-VZ(E)').^2;
    % h2>0 gives a 4x4 Boolean mask eliminating the zero diagonal entries
    % h2(h2>0) is then a (column) vector of the non-zero squared distances
    h2 = min(h2(h2>0));
    if t==1 || h2 < h2_min
        h2_min = h2;
    end
end
param.h_min = sqrt(h2_min);

% Mow call problem config stage 2, to set up elt_to_domain mapping elements
% to domain codes, and phys() structs containing physical properties
% of domains.

param = problem(param.POST_MESH_LOAD, param);

% Throw an error if any of these fields hasn't been assigned in param

check_mandatory(param, {'elt_to_domain', 'phys'});

% Depending on its source, param.elt_to_domain might be a column instead of 
% a row, in which case transpose it. While we're here, check the length
% is correct.

if ~isvector(param.elt_to_domain)
    error('param.elt_to_domain must be one-dimensional');
elseif numel(param.elt_to_domain) ~= K
    error('param.elt_to_domain length %d but %d elements', ...
        numel(param.elt_to_domain), K);
elseif iscolumn(param.elt_to_domain)
    param.elt_to_domain = param.elt_to_domain';
end

% Also check that param.phys is the right shape ...

if ~isvector(param.phys) || ~isstruct(param.phys) || isempty(param.phys)
    error('param.phys must be a non-empty one-dimensional struct array');
end

% ... and has valid domain codes

for domain = unique(param.elt_to_domain)
    if domain ~= round(domain) || domain < 1 || domain > numel(param.phys)
        error('Domain code %d in param.elt_to_domain not a valid index in param.phys');
    end
end

% insert into param struct the locations of elements of different model 
% types, interface faces etc.

param = setup_domains(param);

% Now we know if there are any poroelastic elements, we can decide to use a
% pure elastic or mixed elastic-poroelastic model (unless it's already 
% been set in the problem file, in which case it's non-empty and should be 
% left alone)

if isempty(param.pure_elastic)
    param.pure_elastic = isempty(param.poroelastic_elt);
end

% Identify numbers and storage locations of fields, which differ in the 
% elastic and poroelastic models. In a mixed model, number of elastic 
% fields = 9 is inflated to 13 in DG solver by inserting zeros into \zeta 
% and vfx, vfy, vfz. The matrices Q_e, A_e, B_e, C_e, Qi_e etc. are
% similarly inflated.

% These are the same in both 9-field elastic and 13-field poroelastic models
% First column is for non-adjoint, velocity-strain formulation; second is for
% adjoint, velocity-stress formulation

param.fld_e11  = 1; param.fld_s11 = 1; 
param.fld_e22  = 2; param.fld_s22 = 2; 
param.fld_e33  = 3; param.fld_s33 = 3;
param.fld_e12  = 4; param.fld_s12 = 4;
param.fld_e23  = 5; param.fld_s23 = 5;
param.fld_e13  = 6; param.fld_s13 = 6;

% These differ in 9-field elastic and 13-field poroelastic models

if param.pure_elastic
	param.fld_vx = 7;
	param.fld_vy = 8;
	param.fld_vz = 9;
	param.elastic_fields = 1:9;
	param.num_str_fields = 6;
	param.num_vel_fields = 3;
    param.Nfields = 9;
else
	param.fld_zeta = 7; param.fld_p_f = 7;
	param.fld_vx = 8;
	param.fld_vy = 9;
	param.fld_vz = 10;
	param.fld_vfx  = 11;
	param.fld_vfy  = 12;
	param.fld_vfz  = 13;
	param.elastic_fields = [1:6, 8:10];
	param.num_str_fields = 7;
	param.num_vel_fields = 6;
    param.Nfields = 13;
end

% field locations and storage for sensitivity kernels

if param.run_adjoint_method

    param.num_kernels = 7;

    param.fld_ker_rho_a    = 1;
    param.fld_ker_rho_f    = 2;
    param.fld_ker_m        = 3;
    param.fld_ker_kappa_fr = 4;
    param.fld_ker_mu_fr    = 5;
    param.fld_ker_alpha    = 6;
    param.fld_ker_M        = 7;

    param.sens_ker = zeros(Np, K, param.num_kernels);
    
end

% Copy physical parameters to all elements in each domain and assemble
% mass matrix, Jacobian matrices, compute wavespeeds and associated 
% eigenvectors and the adjoint formulation counterparts

param = build_param(param);
if param.run_forward_adj
    param = build_param_adj(param);
end

% At this point, we expect all the boundary conditions mapX and vmapX
% vectors to be empty. This will not be true if BCType has non-zero entries
% (most likely, it's empty at this point, which is fine)
% TODO this should be updated so it is possible to read BCType from mesh 
% generator output; see comments above, where mesh data is read and below,
% where BCType is reverse engineered from mapX data

if any(BCType(:) ~= 0)
    error('Unexpected non-zero values in BCType');
end

% Now call problem config stage 3, to set up boundary conditions, sources,
% receivers and initial values. This needs to happen after domain maps 
% have been set up and physical parameters have been mapped to elements

param = problem(param.POST_DOMAIN_SETUP, param);

% Check consistency of vmapX with mapX

if any(vmapC ~= vmapM(mapC))
    error('vmapC not consistent with mapC')
elseif any(vmapD ~= vmapM(mapD))
    error('vmapD not consistent with mapD')
elseif any(vmapF ~= vmapM(mapF))
    error('vmapF not consistent with mapF')
elseif any(vmapI ~= vmapM(mapI))
    error('vmapI not consistent with mapI')
elseif any(vmapN ~= vmapM(mapN))
    error('vmapN not consistent with mapN')
elseif any(vmapO ~= vmapM(mapO))
    error('vmapO not consistent with mapO')
elseif any(vmapS ~= vmapM(mapS))
    error('vmapS not consistent with mapS')
elseif any(vmapW ~= vmapM(mapW))
    error('vmapW not consistent with mapW')
end

% Create BCType from boundary maps above. This is not the usual
% way round, but it allows us to save BCType to be read later
% by the C framework. 
% TODO this should be updated so it is possible to read BCType from mesh 
% generator output; see comments above, where mesh data is read, and
% where an error is thrown if BCType has non-zero entries before problem
% config stage 3 is called.

BCType = zeros(Nfaces*K, 1);

BCType(deflate(mapC, Nfp)) = Cyl;
BCType(deflate(mapD, Nfp)) = Dirichlet;
BCType(deflate(mapF, Nfp)) = Far;
BCType(deflate(mapI, Nfp)) = In;
BCType(deflate(mapN, Nfp)) = Neuman; % spelling from H&W
BCType(deflate(mapO, Nfp)) = Out;
BCType(deflate(mapS, Nfp)) = Slip;
BCType(deflate(mapW, Nfp)) = Wall;

BCType = reshape(BCType, [Nfaces, K])';

% from param.g_rec, the coordinates of the receivers, find the elements
% and interpolating operators

param = setup_receiver(param);
param.num_rec = size(param.g_rec,1);

% In the adjoint method problem, these receivers becomes sources;
% set these up if we are running this.

if param.run_adjoint_method
    param = setup_adjoint_wavefield_source(param);
end

% Time stepping. Need maximal wave speed for CFL condition
% Note cp_ela and cpI_i are treated differently in build_param:
% cp_ela has only elastic elements, cpI_i has all elements, hence
% different-looking code. TODO Fix this at some point.
% Absolute values required because wavespeeds are negative if
% integrating backwards in time

if ~isempty(param.elastic_elt)
    v_max = max(abs(param.cp_ela));
else
    v_max = max(abs(param.cpI_i(param.poroelastic_elt)));
end
fprintf("Maximal wave speed %.5g\n", v_max);
fprintf("Source frequency %.5g\n", param.freq_sou);

% CFL condition. h_min was calculated when mesh was read, param.cfl
% defaults to 0.4 but can be overridden in problem config file

param.dt = param.cfl*param.h_min/(N^2*v_max);

% save data for C system, if required

if param.save_params

    % Reshape the initial values (param.field and, if there are HF
    % elements, param.mem_hf) to the form meeded by the C system. Save 
    % temporarily as param.initial; this will be removed from
    % param when the data has been saved.

    param.initial = combine_memory(param.field, param.mem_hf, mem_hf, 2, param, true);

    % If we have an exact solution, evaluate it at the final time and
    % save temporarily as param.exact_final; this will be removed from
    % param when the data has been saved. Note that this uses
    % initial_time and time_width as specified in the problem config file;
    % it does not take account of save/restart or manual setting of 
    % param.num_steps
    
    if isfield(param, 'exact')
       
        if isempty(param.hf_elt)
            field_exact = param.exact(param.initial_time+param.time_width, x, y, z, param);
            param.exact_final = reshape(field_exact, K*Np, param.Nfields);
        else
            [field_exact, mem_hf_exact] = ...
                param.exact(param.initial_time+param.time_width, x, y, z, param);
            param.exact_final = combine_memory(field_exact, mem_hf_exact, 2, param, true);
        end
    end

    % Now save the whole param structure and the reverse-engineered 
    % BCType datafor the C system to find.
    
    struct_to_h5(param, fullfile(param.problem_data_dir, param.C_extra_filename)); 
    dlmwrite(fullfile(param.problem_data_dir, bctype_filename), BCType, ' ');
    
    % Remove some unnecessary clutter from the param struct once it's been 
    % saved: the Matlab code can use the function handle param.exact to
    % find exact soution at any time, and param.initial is just a 
    % repackaging of param.field and param.mem_hf, which are
    % used by the Matlab system
    
    param = rmfield(param, 'initial');
    if isfield(param, 'exact_final')
        param = rmfield(param, 'exact_final');
    end
    
end % if param.save_params % (save for C system)

% If required, load saved state from earlier run and use it as the initial
% state for this run.
% Note that the time loop always runs from time 0 and step 0 but 
% offset_time and offset_step are added when time is used in e.g. 
% source terms or boundary conditions and step is used in e.g receiver 
% storage. After each run, the final time (including offset) is saved, to 
% be loaded as the offset for the next run. frame is the index number of 
% the run within the current sequence, used for saving a sequence of 
% visualisations

if param.restart_from_saved
    h5_file = fullfile(param.problem_data_dir, param.save_state_filename);

    % times of each RK step in run to this point
    
    saved_step_time = h5read(h5_file, '/step_time');
    
    % final time of saved simulation becomes initial time of new
    % simulation
    
    param.initial_time = saved_step_time(end);
    
    % field and mem_hf are the last state of the saved run; receivers are
    % all values recorded up to the end of the saved run. 
    % h5load is a custom function that returns empty for non-existent
    % datasets, avoided the problem that h5create won't create them

    param.field        = h5load(h5_file, '/field');
    param.field_adj    = h5load(h5_file, '/field_adj');
    param.mem_hf       = h5load(h5_file, '/mem_hf');
    param.mem_hf_adj   = h5load(h5_file, '/mem_hf_adj');
    param.sens_ker     = h5load(h5_file, '/sens_ker');
    saved_receiver     = h5load(h5_file, '/receiver');
    saved_receiver_adj = h5load(h5_file, '/receiver_adj');
    
    % count number of times simulation has been (re)started, used for 
    % e.g. naming figure files
    
    frame = h5read(h5_file, '/frame') + 1;

else % Not restarting from saved
    
    frame = 0;
    saved_receiver = [];
    saved_receiver_adj = [];
    saved_step_time = [];

end % param.restart_from_saved

% Finally, we can run the numerics!

if param.run_simulation
    
    problem(param.PRE_SIMULATION, param);
    
    % Run time-stepping loop
    output = biot_3d(param);

    problem(param.POST_SIMULATION, param, output);

    % used for plotting and comparison with exact solution
    final_time = output.step_time(end);

    if param.pure_elastic
        fname = {'e11 ', 'e22 ', 'e33 ', 'e12 ', 'e23 ', 'e13 ', ...
            'vx  ', 'vy  ', 'vz  '};
    else
        fname = {'e11 ', 'e22 ', 'e33 ', 'e12 ', 'e23 ', 'e13 ', 'zeta', ...
            'vx  ', 'vy  ', 'vz  ', 'vfx ', 'vfy ', 'vfz '};
    end
    
    if isfield(param, 'exact')
        % If there are HF elements, param.exact returns the memory
        % variables as separate outputs; calling with exact_all as the only
        % destination variable effectively discards the memory variable
        % output.
        exact_all = param.exact(final_time, x, y, z, param);
        fprintf('Field L^2 exact     L^2 simul     L^2 error     L^2 Relative error\n');
        for fld = 1:param.Nfields
            simul = output.field(:,:,fld);
            exact = exact_all(:,:,fld);
            fprintf("%-5s %.7e %.7e %.7e %.7e\n", fname{fld}, L2norm(exact), ...
            L2norm(simul), L2norm(exact-simul), L2norm(exact-simul)/L2norm(exact));
        end
    end % if isfield(param, 'exact')
end  % if param.run_simulation

% if we have loaded a previously saved array of receiver values, prepend
% these to the values from the run (receiver arrays are 3-dimensional, 
% indexed by (receiver, step, field). Note that receiver values are also
% saved for the initial data, so the last values of the saved data should
% match the first values of the new data, and one of these duplicates
% should be discarded.
% TODO Does this behave properly if run_simulation is false and
% restart_from_saved is true? It should replicate the saved state at this
% point.

if ~isempty(saved_step_time)
    if saved_step_time(end) ~= output.step_time(1)
        fprintf('Saved step times values do not match newly created\n')
    end
    output.step_time = [saved_step_time, output.step_time(2:end)];
end
if ~isempty(saved_receiver)
    if ~isequal(saved_receiver(:,end,:), output.receiver(:,1,:))
        fprintf('Saved receiver values do not match newly created\n')
    end
    output.receiver = cat(2, saved_receiver, output.receiver(:,2:end,:));
end
if ~isempty(saved_receiver_adj)
    if ~isequal(saved_receiver_adj(:,end,:), output.receiver_adj(:,1,:))
        fprintf('Saved receiver_adj values do not match newly created\n')
    end
    output.receiver_adj = cat(2, saved_receiver_adj, output.receiver_adj(:,2:end,:));
end

% Save final state for later restart, if required

if param.save_for_restart
    
    % Need absolute path here, because exist searches the matlab path
    
    h5_file = fullfile(param.problem_data_dir, param.save_state_filename);
    if exist(h5_file, 'file') == 2
        delete(h5_file);
    end

    % h5save is a custom routine that calls h5create and h5write; it also
    % silently skips empty arrays, which h5create won't handle
    
    h5save(h5_file, '/step_time',    output.step_time);
    h5save(h5_file, '/field',        output.field);
    h5save(h5_file, '/field_adj',    output.field_adj);
    h5save(h5_file, '/mem_hf',       output.mem_hf);
    h5save(h5_file, '/mem_hf_adj',   output.mem_hf_adj);
    h5save(h5_file, '/receiver',     output.receiver);
    h5save(h5_file, '/receiver_adj', output.receiver_adj);
    h5save(h5_file, '/sens_ker',     output.sens_ker);

    % Counts the number of times the simulation has been (re)started
    
    h5save(h5_file, '/frame', frame);
    
end

% Save output for comparison with C output, if required.

if param.save_output_field || param.save_output_rec
    % Need absolute path because exist searches whole path
    h5_file = fullfile(param.problem_data_dir, param.output_data_filename);
    if exist(h5_file, 'file') == 2
        delete(h5_file);
    end
end

if param.save_output_field
    
    % Combine main fields with memory variables (if they exist).
    % These will be empty if the corresponding output fields are empty.

    save_field     = combine_memory(output.field,     output.mem_hf,     2, param, true);
    save_field_adj = combine_memory(output.field_adj, output.mem_hf_adj, 2, param, true);

    h5save(h5_file, '/field', save_field);
    h5save(h5_file, '/field_adj', save_field_adj);
    
end % if param.save_output_field

if param.save_output_rec
    
    h5save(h5_file, '/receiver', output.receiver);
    h5save(h5_file, '/receiver_adj', output.receiver_adj);

end % if param.save_output_rec

% If param.save_adjoint_data is set, save the receiver values in the file
% specified by param.save_adjoint_filename. This has three datasets, 
% '/receiver_data' (measurement data), '/receiver_adj_model' (adjoint model 
% receiver values) and 'field_adj_model' (adjoint model final wavefield). 
% measurement data must be saved first; if the h5 file exists, it is 
% deleted and recreated empty, then '/receiver_data' is created and 
% populated. '/receiver_adj_model' and '/field_adj_model' are created 
% (if they do not exist) and populated on subsequent runs and can be 
% overwritten any number of times, as long as the dataset sizes do not
% change.

if param.save_adjoint_data == param.SAVE_RECEIVER_DATA
    % Need absolute path because exist searches whole path
    h5_file = fullfile(param.problem_data_dir, param.save_adjoint_filename);
    if exist(h5_file, 'file') == 2
        delete(h5_file);
    end
    h5create(h5_file, '/receiver_data', size(output.receiver), 'datatype', 'double');
    h5write(h5_file, '/receiver_data', output.receiver);
    h5create(h5_file, '/step_time', size(output.step_time), 'datatype', 'double');
    h5write(h5_file, '/step_time', output.step_time);
end
    
if param.save_adjoint_data == param.SAVE_ADJ_MODEL
    h5_file = fullfile(param.problem_data_dir, param.save_adjoint_filename);
    try
        h5create(h5_file, '/receiver_adj_model', size(output.receiver), 'datatype', 'double');
        h5create(h5_file, '/field_adj_model', size(output.field), 'datatype', 'double');
    catch
    end
    h5write(h5_file, '/receiver_adj_model', output.receiver);
    h5write(h5_file, '/field_adj_model', output.field);
end

end % function driver_biot_3d(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mandatory parameters without defaults: keynames should be a 1-dimensional
% cell array. If any of these fields have not been assigned in the struct
% param, throw an error. Intended to pick up problems at the earliest 
% possible point, to give a more intelligible error message. Note use of 
% keyname{1}: keyname itself is not a string, but a 1x1 cell array 
% containing a string.

function check_mandatory(param, keynames)
    for keyname = keynames
        if ~isfield(param, keyname{1})
            error('Problem file has not specified param.%s', keyname{1})
        end
    end
end
