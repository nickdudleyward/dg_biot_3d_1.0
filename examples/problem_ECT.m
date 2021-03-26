function param = problem_ECT(setup_stage, param, output)

% Problem file for Elastic Convergence Test
% SPE 2019/11/29, updated 2020/07/21

% Hesthaven-Warburton global variables
Globals3D;

% Setup happens in three separate stages. This routine is called three
% times, with different values of setup_stage.

switch setup_stage
    
    case param.PRE_MESH_LOAD
        % Here we must specify the mesh and polynomial order. Lots of other
        % global parameters can also go here, but they all have defaults;
        % of these, we only set param.exact to a handle to an exact solver
        % (which is used only for status messages in the Matlab code).
        param.mesh = 'meshes/mesh_0';
        param.order = 4;
        param.exact = @plane_wave_ie;
        param.problem_data_dir = '~/DG_run_3d/ECT'; % change as appropriate
        param.pure_elastic = false;
        %param.num_steps = 10;
    
    case param.POST_MESH_LOAD
        % Here we must assign domain codes to elements (unless they were
        % read from a domain index file) and physical parameters
        % to domain codes. param.phys() could actually be set up in the
        % previous config stage, but elt_to_domain must go here, and it
        % seems sensible to keep elt_to_domaim and phys() together.
        
        % All elements have domain code 1 ...
        param.elt_to_domain = ones(1, K);
        % ... which has an elastic model type and elastic parameter values
        % Build domain 1 setup in local struct p ...
        p.model_type = param.ELASTIC;
        p.rho_ela = 3500;
        p.cp_ela = 6000;
        p.cs_ela = 3000;
        p.mu_ela = p.cs_ela^2*p.rho_ela;
        p.lambda_ela = p.cp_ela^2*p.rho_ela-2*p.mu_ela;
        % ... then assign to param.phys(1), thus avoiding tedious 
        % repetition of param.phys(1). above
        param.phys(1) = p;
    
    case param.POST_DOMAIN_SETUP
        % Here we set up boundary conditions, sources, receivers and initial 
        % values. In this case, much of this could have happened in the
        % previous stage, but initial values have to go here; in more
        % complicated examples, source setup has to go here (?). 
        % BC setup goes here in case BC type of boundary face depends 
        % on model type of boundary element (?)
        
        % The whole boundary is subject to Dirichlet conditions
        mapD  = mapB;
        vmapD = vmapB;
        % boundary value function handle for the Matlab code; CT_model_type
        % is the equivalent setting for the C code
        param.DirichletBC = @plane_wave_ie;
        param.CT_model_type = param.ELASTIC;
        % no sources (or receivers; this is the default state)
        % (source specification might change in the near future!)
        param.elm_sou = [];
        % Initial values. No need to set param.[ex,ey,ez] as these
        % default to empty, which is correct when there are no HF elements
        param.field = plane_wave_ie(0, x, y, z, param);

    case { param.PRE_SIMULATION, param.POST_SIMULATION }
        % No action in this problem file
    
    otherwise
        error('Invalid setup stage code %d', setup_stage);
        
end
