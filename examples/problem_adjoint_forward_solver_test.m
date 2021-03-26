function param = problem_adjoint_forward_solver_test(setup_stage, param, output)

% Problem file for forward adjoint problem
% SPE/NFDW 2020/07/07, updated 2020/07/21

% Hesthaven-Warburton global variables
Globals3D;

% Setup happens in three separate stages. This routine is called three
% times, with different values of setup_stage.

switch setup_stage
    
    case param.PRE_MESH_LOAD
        % Here we must specify the mesh and polynomial order. Lots of other
        % global parameters can also go here, but they all have defaults.
        param.mesh = 'meshes/mesh_0';
        param.order = 2;
        param.problem_data_dir = '~/DG_run_3d/forward_solver_test'; % change as appropriate
        % These must be set before build_param is called; might as well
        % be here
        param.run_forward = true;
        param.run_forward_adj = true;
        param.num_steps = 100; % RK steps; overrides any param.time_width setting
        %param.time_width = 1/param.freq_sou;
        param.steps_per_report = 10;
    
    case param.POST_MESH_LOAD

        % Here we must assign domain codes to elements (unless they were
        % read from a domain index file) and physical parameters
        % to domain codes. param.phys() could actually be set up in the
        % previous config stage, but elt_to_domain must go here, and it
        % seems sensible to keep elt_to_domaim and phys() together.

        % Elements below this z value will form domain 1 (UPPER INVISCID). 
        % Elements above will form domain 2 (LOWER INVISCID)
        z0 = 2.3;
        % To identify excat boundary in a regular mesh, snap to nearest 
        % z coordinate of element vertex (not compuational node)
        [~, j] = min(abs(z0-VZ));
        z0 = VZ(j);
        fprintf('Snapped to z=%g\n', z0);
        % To assign model type by location, find centroid of each element 
        %x_mean = mean(VX(EToV),2);
        %y_mean = mean(VY(EToV),2);
        z_mean = mean(VZ(EToV),2);
        % Can now assign domains based on element location:
        % domain 1 has centroid below z0, domain 2 above
        param.elt_to_domain(z_mean <= z0)  = 1; % model type INVISCID
        param.elt_to_domain(z_mean >  z0)  = 2; % model type INVISCID
                
        % Build domain 1 (UPPER INVISCID) setup in local struct p1
        p1.model_type = param.INVISCID;
        p1.rho_s  = 2200;
        p1.rho_f  = 950;
        % poroelastic moduli
        p1.mu_fr = 3e9;
        p1.k_s   = 6.9e9;
        p1.k_f   = 2e9;
        p1.k_fr  = 6.7e9;
        % fluid and porous media
        p1.eta = 0.0; %viscosity
        p1.phi = 0.4;%porosity
        p1.tau = 2;%tortuosity

        % Build domain 2 (LOWER INVISCID) setup in local struct p2
        % densities
        p2.model_type = param.INVISCID;
        p2.rho_s  = 2650;
        p2.rho_f  = 750;
        % poroelastic moduli
        p2.mu_fr = 4.4e9;
        p2.k_s   = 37e9;
        p2.k_f   = 1.7e9;
        p2.k_fr  = 2.2e9;
        % fluid and porous media
        p2.eta = 0.0; %viscosity
        p2.phi = 0.2;%porosity
        p2.tau = 2;%tortuosity
        
        % Assign to param.phys(1,2) (avoiding tedious repetition of 
        % param.phys(1). above). Note use of custom function concat_struct
        % which allow concatenation of struct arrays with different fields.
        
        param.phys = concat_struct(p1, p2);
            
    case param.POST_DOMAIN_SETUP
        % Here we set up boundary conditions, sources, receivers and initial 
        % values. In this case, much of this could have happened in the
        % previous stage, but initial values have to go here; in more
        % complicated examples, source setup has to go here (?). 
        % BC setup goes here in case BC type of boundary face depends 
        % on model type of boundary element (?)
        
        % for a cuboid mesh, set the top face BCs (i.e. those with outward
        % normal (0,0,1) to FP precision) to W=Wall and all other
        % BCs to outflow=Far=F

        onz = nz(:);
        tol = 1e-10;

        mapW  = intersect(mapB, find(abs(onz-1)<tol));
        vmapW = vmapM(mapW);
        mapF  = setdiff(mapB, mapW);
        vmapF = vmapM(mapF);
        
        % Tell the C code that there is no Dirichlet BC function; 
        % for Matlab, this corresponds with leaving param.DirichletBC unset
        param.CT_model_type = param.NONE;
        % Receiver location(s) are rows of param.g_rec
        param.g_rec = [0.75 0.75 2];
 
        % Set element containing source, r, s,t coordinates of source
        % within element, model type of element within param
        % (source specification might change in the near future!)
        
        %param.source_type = param.SURFACE_FORCE;
        param.source_type = param.MOMENT_TENSOR;
        
        switch param.source_type
            case param.MOMENT_TENSOR
                % Source location, moment tensor, Ricker wavelet frequency        
                param.g_sou = [2.5, 2.5, 3.5];
                param.t0_sou = 1.2/param.freq_sou;
                param.M0 = 1e6;
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
        
        
        % Initial values: identically zero. No need to set param.[ex,ey,ez] 
        % as these default to empty, which is correct when there are no 
        % HF elements
        param.field = zeros(Np, K, param.Nfields);
        param.field_adj = zeros(Np, K, param.Nfields);

    case { param.PRE_SIMULATION, param.POST_SIMULATION }
        % No action in this problem file
    
    otherwise
        error('Invalid setup stage code %d', setup_stage);
        
end
