function MGC = setup_magic()

% These are magic numbers used in various parts of the DG system for Biot
% equations. There should be a clean way to keep these synchronised
% between Matlab, C and Python, but for the time being this is just for
% Matlab

MGC = struct();

% Source types

MGC.MOMENT_TENSOR = 1;
MGC.SURFACE_FORCE = 2;

% setup stages

MGC.PRE_MESH_LOAD     = 1;
MGC.POST_MESH_LOAD    = 2;
MGC.POST_DOMAIN_SETUP = 3;
MGC.PRE_SIMULATION    = 4;
MGC.POST_SIMULATION   = 5;

% model types. first digit is main type (elastic / poro), second is 
% subtype (for poro, inviscid, low freq, high freq)
% NB param.POROELASTIC is a dummy type code, used only to distinguish 
% elastic types ( < param.POROELASTIC) and poroelastic types
% ( > param.POROELASTIC).
% code NONE = 0 is used by the C framework; applied to the whole mesh,
% it means that it's not all of one model type

MGC.NONE = 0;
MGC.ELASTIC = 11;
MGC.POROELASTIC = 20;
MGC.INVISCID = 21;
MGC.LOW_FREQUENCY = 22;
MGC.HIGH_FREQUENCY = 23;

% Operator splitting methods

MGC.SPL_NONE = 0;
MGC.SPL_SIMPLE_AN = 1; 
MGC.SPL_SIMPLE_NA = 2; 
MGC.SPL_STRANG_ANA = 3;

% Runge-Kutta methods

MGC.RK_LSEX = 1;
MGC.RK_LSIMEX = 2;

% Saving receiver values for adjoint method

MGC.SAVE_ADJOINT_NONE = 0;
MGC.SAVE_RECEIVER_DATA = 1;
MGC.SAVE_ADJ_MODEL = 2;

