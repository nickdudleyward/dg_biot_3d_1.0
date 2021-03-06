function param = problem_adjoint_method_test_2(setup_stage, param, output)

% Do not run this file directly; run driver_test_new with this file's name
% passed as a parameter or assigned to variable default_problem_file at the
% top of driver_test_new

% Problem file for adjoint method problem solved in paper

% SPE/NFDW 2020/07/14, updated 2020/07/21

% Magic constants for later in the run, distinguishing three different
% modes of operation. Set run_type to one of these in case 
% param.PRE_MESH_LOAD; this will be remembered for subsequent cases
% because run_type is declared persistent.

persistent run_type;
MAKE_ADJOINT_DATA = 1;
SETUP_ADJOINT_PROBLEM = 2;
SOLVE_ADJOINT_METHOD = 3;

% Hesthaven-Warburton global variables

Globals3D;

% Setup usually happens in three separate stages; this routine is called 
% once per stage, with different values of setup_stage. There is also a
% post-simulation stage, for output processing.

switch setup_stage
    
    case param.PRE_MESH_LOAD
        
        % enable one of these only!
        %run_type = MAKE_ADJOINT_DATA;
        %run_type = SETUP_ADJOINT_PROBLEM;
        run_type = SOLVE_ADJOINT_METHOD;
        
        % Here we must specify the mesh and polynomial order. Lots of other
        % global parameters can also go here, but they all have defaults.
        param.mesh = 'meshes/mesh_0';
        param.order = 4;
        
        param.problem_data_dir = '~/DG_run_3d/adjoint_method_test_2'; % change as appropriate
        param.save_adjoint_filename = 'biot_3d_adjoint_rec_data.h5';

        % default; can be overridden in switch below
        param.time_width = 6/param.freq_sou;
        
        switch run_type
            
            case MAKE_ADJOINT_DATA
                param.run_forward = true;
                param.run_forward_adj = false;
                param.run_adjoint_method = false;
                %param.num_steps = 10; % #RK steps; overrides any param.time_width setting
                param.time_width = 6/param.freq_sou; % make explicit 6
                param.steps_per_report = 5;
                param.save_adjoint_data = param.SAVE_RECEIVER_DATA;
        
            % Generate final field and receiver values to assemble
            % adjoint wavefield source using adjoint data
            case SETUP_ADJOINT_PROBLEM
                param.restart_from_saved = false;
                param.save_for_restart = false;
                param.run_forward = true;
                param.run_forward_adj = false;
                param.run_adjoint_method = false;
                param.save_adjoint_data = param.SAVE_ADJ_MODEL;
                %param.num_steps = 1; % RK steps; overrides any param.time_width setting
                param.time_width = 6/param.freq_sou; % 6
                param.steps_per_report = 5;

            case SOLVE_ADJOINT_METHOD
                % param.initial_time is set in setup_adjoint_wavefield_source
                % (and cannot be overridden here as that is called after
                % POST_DOMAIN_SETUP!)
                param.restart_from_saved = false;
                param.save_for_restart = true;
                param.run_forward = true;
                param.run_forward_adj = true;
                param.run_adjoint_method = true;
                %param.num_steps = -5; % RK steps; overrides any param.time_width setting
                param.time_width = -1/param.freq_sou;
                param.steps_per_report = 5;
                
            otherwise
                error('Invalid run_type %d', run_type)
                
        end

    case param.POST_MESH_LOAD

        % Here we must assign domain codes to elements (unless they were
        % read from a domain index file) and physical parameters
        % to domain codes. param.phys() could actually be set up in the
        % previous config stage, but elt_to_domain must go here, and it
        % seems sensible to keep elt_to_domain and phys() together.

        % The following parameters are used to generate the receiver data
        % for the adjoint method example and will be used to build the 
        % adjoint source
        
        % Select one element to have different parameter values 
        % corresponding to p2 = phys(2) below
        param.g_jumps = [2.5, 2.5, 4];
        [param.elm_jumps, ~] = locate_element(param.g_jumps);
        
        % Build domain 1 (all but one element) setup in local struct p1
        p1.model_type = param.INVISCID;
        p1.rho_s  = 2650;
        p1.rho_f  = 900;
        % poroelastic moduli
        p1.mu_fr = 5e9;
        p1.k_s   = 12e9;
        p1.k_f   = 2e9;
        p1.k_fr  = 10e9;
        % fluid and porous media
        p1.eta = 0.0; %viscosity
        p1.phi = 0.3;%porosity
        p1.tau = 1.2;%tortuosity

        % Build domain 2 (one exceptional element) setup in local struct p2
        % densities
        p2 = p1;
        p2.rho_s = 2*p2.rho_s;
        p2.rho_f = 2*p2.rho_f;
 
        % Depending on run_type, assign to param.phys(1,2) (avoiding 
        % tedious repetition of physical parameters above). Note use of 
        % custom function concat_struct which allows concatenation of 
        % struct arrays with different fields.
        
        switch run_type
            case MAKE_ADJOINT_DATA
                param.elt_to_domain = ones(1, K);
                param.elt_to_domain(param.elm_jumps) = 2;
                param.phys = concat_struct(p1, p2);
            case { SETUP_ADJOINT_PROBLEM, SOLVE_ADJOINT_METHOD }
                param.elt_to_domain = ones(1, K);
                param.phys = p1;
            otherwise
                error('Invalid run_type %d', run_type)            
        end
        
    case param.POST_DOMAIN_SETUP
        % Here we set up boundary conditions, sources, receivers and initial 
        % values. In this case, much of this could have happened in the
        % previous stage, but initial values have to go here; in more
        % complicated examples, source setup has to go here (?). 
        % BC setup goes here in case BC type of boundary face depends 
        % on model type of boundary element (?)
        
        % Whole boundary is reflecting for adjoint problem

        mapW  = mapB;
        vmapW = vmapM(mapW);
        
        % Tell the C code that there is no Dirichlet BC function; 
        % for Matlab, this corresponds with leaving param.DirichletBC unset
        param.CT_model_type = param.NONE;
        
        % Receiver location(s) are rows of param.g_rec
        
        num_rec = 100;
        rec_line = linspace(0,5,sqrt(num_rec));

        [X,Y] = meshgrid(rec_line);
        X = X(:);
        Y = Y(:);
        param.g_rec = [X Y 5*ones(num_rec,1)];
        
        % Set element containing source, r, s,t coordinates of source
        % within element, model type of element within param
        % (source specification might change in the near future!)
        
        %param.source_type = param.SURFACE_FORCE;
        param.source_type = param.MOMENT_TENSOR;
        
        switch param.source_type
            case param.MOMENT_TENSOR
                % Source location, moment tensor, Ricker wavelet frequency        
                param.g_sou = [2.5, 2.5, 2];
                param.t0_sou = 1.2/param.freq_sou;
                param.M0 = 1e2;
                param.Mxx = 1;
                param.Myy = 1;
                param.Mzz = 1;
                param.Mxy = 0;
                param.Mxz = 0;
                param.Myz = 0;
                
                param = locate_source(param);
                param.sou = setup_source_moment(param);
                                  
            case param.SURFACE_FORCE
                param.g_sou = [2.5, 2.5, 5]; %location on surface
                param.t0_sou = 1.2/param.freq_sou;
                param.force = [0, 0, -1e4];
                
                param = locate_source(param);
                param.sou = setup_source_force(param);
                
            otherwise
                error('Invalid source type %d', param.source_type)
        end

        % Setup initial values
        
        switch run_type
            
            case MAKE_ADJOINT_DATA
                param.field = zeros(Np, K, param.Nfields);
                
            case SETUP_ADJOINT_PROBLEM
                param.field = zeros(Np, K, param.Nfields);
                
            case SOLVE_ADJOINT_METHOD
                h5_file = fullfile(param.problem_data_dir, param.save_adjoint_filename);
                param.field = h5read(h5_file, '/field_adj_model');
                param.field_adj = zeros(Np, K, param.Nfields);
                
            otherwise
                error('Invalid run_type %d', run_type)            
        end

    case param.PRE_SIMULATION
        % No action in this problem file

    case param.POST_SIMULATION
        
        final_time = output.step_time(end);

        % BEGIN PLOTTING

        switch run_type 
            
            case MAKE_ADJOINT_DATA
                figure;
                set(gcf,'Position',[10 10 800 500]);
                row = 4;
                for i = 1:param.num_rec/10
                    subplot(2,ceil(param.num_rec/20),i)
                    % row*10+i is a receiver number, : is the time sequence
                    % and 8 is field vx
                    plot(output.step_time, output.receiver(row*10+i,:,8))
                    xlim([output.step_time(1) output.step_time(end)])
                    set(gca,'fontsize', 10)
                    title('Measured receiver vx')
                end
                
            case SETUP_ADJOINT_PROBLEM
                figure;
                set(gcf,'Position',[10 10 800 500]);
                row = 4;
                h5_file = fullfile(param.problem_data_dir, param.save_adjoint_filename);
                tmp = output.receiver - h5read(h5_file, '/receiver_data');
                for i = 1:param.num_rec/10
                    subplot(2,ceil(param.num_rec/20),i)
                    % row*10+i is a receiver number, : is the time sequence
                    % and 8 is field vx
                    plot(output.step_time, tmp(row*10+i,:,8))
                    xlim([output.step_time(1) output.step_time(end)])
                    set(gca,'fontsize', 10)
                    title('Adjoint wavefield source vx-d_x')
                end
                
            case SOLVE_ADJOINT_METHOD

                % figure(4)
                % pngfile = sprintf('figs/image_%03d_final_value_pl_order_3.png', frame);
                % print(pngfile, '-dpng', '-r300');

                % extract solid velocity components and compute speed
                f = sqrt(output.field(:,:,8).^2+output.field(:,:,9).^2+output.field(:,:,10).^2);            
                plot_section_3d(2.5, [0 5], [0 5], f, 150);
                title(sprintf('time reversed wave field at %d/f0', final_time*param.freq_sou));

                f = sqrt(output.field_adj(:,:,8).^2+output.field_adj(:,:,9).^2+output.field_adj(:,:,10).^2);            
                plot_section_3d(2.5, [0 5], [0 5], f, 150);
                title(sprintf('time reversed adjoint wave field at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_rho_a), 150);
                title(sprintf('ker\\_rho\\_a at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_rho_f), 150);
                title(sprintf('ker\\_rho\\_f at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_m), 150);
                title(sprintf('ker\\_m at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_kappa_fr), 150);
                title(sprintf('ker\\_kappa\\_fr at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_mu_fr), 150);
                title(sprintf('ker\\_mu\\_fr at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_alpha), 150);
                title(sprintf('ker\\_alpha at %d/f0', final_time*param.freq_sou));

                plot_section_3d(2.5, [0 5], [0 5], output.sens_ker(:,:,param.fld_ker_M), 150);
                title(sprintf('ker\\_M at %d/f0', final_time*param.freq_sou));

        end

        % END PLOTTING                
    otherwise
        error('Invalid setup stage code %d', setup_stage);
        
end
