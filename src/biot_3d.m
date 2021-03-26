% Modified for forward adjoint case 2020/07/09 SPE/NFDW
% Combined memory variables into single array 2020/07/30 SPE/NFDW
% Operator splitting, IMEX, dissipation not yet implemented for adjoint

% This codes implements the DG solver for Biot model for the inviscid, low 
% frequency and high frequency regimes.

function output = biot_3d(param)

% copies of input data, to be modified and returned as output
    
field = param.field;
mem_hf = param.mem_hf;
field_adj = param.field_adj;
mem_hf_adj = param.mem_hf_adj;
sens_ker = param.sens_ker;

% sensitivity kernels, awaiting adjoint method implementation

ker_rho_a = [];
ker_kappa_fr = [];

Globals3D;

fprintf('Main loop\n');

% Number of time steps.

% If param.num_steps is non-zero, its absolute value is the number of time
% steps taken and its sign is the direction of integration. The time step
% is param.dt, from the CFL condition.
%
% If param.num_steps is zero, the absolute value of param.time_width is the
% width of the integration region and its sign is the direction of 
% integration. The time step is reduced from param.dt so an integer number
% of equal time steps are used.
%
% param.dt is always positive. dt is positive or negative, according to the
% direction of integration. dt is typically smaller in magnitude than
% param.dt if param.num_steps is not specified, or equal in magnitude if 
% param.num_steps is specified.

initial_time = param.initial_time;
time_width = param.time_width;

if param.num_steps ~= 0
    dt = sign(param.num_steps)*param.dt;
    num_steps = abs(param.num_steps);
    time_width = num_steps*dt;
else
    num_steps = ceil(abs(time_width/param.dt));
    dt = time_width / num_steps;
end

fprintf('Number of time steps %d\n', num_steps);
fprintf('dt: %e\n', dt)

%keyboard
% to check the CPU time and wallclock time
t0 = cputime;
tic;
if param.pure_elastic
    fname = {'e11 ', 'e22 ', 'e33 ', 'e12 ', 'e23 ', 'e13 ', ...
        'vx  ', 'vy  ', 'vz  '};

    fname_adj = {'s11 ', 's22 ', 's33 ', 's12 ', 's23 ', 's13 ', ...
            'vx  ', 'vy  ', 'vz  '};
else
    fname = {'e11 ', 'e22 ', 'e33 ', 'e12 ', 'e23 ', 'e13 ', 'zeta', ...
        'vx  ', 'vy  ', 'vz  ', 'vfx ', 'vfy ', 'vfz '};

    fname_adj = {'s11 ', 's22 ', 's33 ', 's12 ', 's23 ', 's13 ', 'p_f ', ...
            'vx  ', 'vy  ', 'vz  ', 'vfx ', 'vfy ', 'vfz '};
end

% Extract coefficient data from RK.h5 for the method specified by the
% string param.RK_method (which should correspond to a group name in the
% h5 file).
% This file is in the src directory which is on the path, so no
% absolute path needed (but a copy in pwd will take priority)
% If no method is specified, fall back on the method used in Hesthaven 
% and Warburton.

h5_RK = 'RK.h5';

if isfield(param, 'RK_method')
    prefix = strcat('/', param.RK_method, '/');
    RK_type = h5read(h5_RK, strcat(prefix, 'RK_type'));
    if RK_type == param.RK_LSEX
        LSEXa = h5read(h5_RK, strcat(prefix, 'LSEXa'));
        LSEXb = h5read(h5_RK, strcat(prefix, 'LSEXb'));
        LSEXc = h5read(h5_RK, strcat(prefix, 'LSEXc'));
        num_stages = numel(LSEXa);
    elseif RK_type == param.RK_LSIMEX
        LSIMEXaIMd  = h5read(h5_RK, strcat(prefix, 'LSIMEXaIMd'));
        LSIMEXaIMsd = h5read(h5_RK, strcat(prefix, 'LSIMEXaIMsd'));
        LSIMEXaEXsd = h5read(h5_RK, strcat(prefix, 'LSIMEXaEXsd'));
        LSIMEXb     = h5read(h5_RK, strcat(prefix, 'LSIMEXb'));
        LSIMEXc     = h5read(h5_RK, strcat(prefix, 'LSIMEXc'));
        num_stages  = numel(LSIMEXaIMd);
    else
        error('Unknown RK_type code %d', RK_type)
    end
else % param.RK_method not set; fall back on 5-stage order 4 LSRK default
    % from Hesthaven-Warburton library globals
    RK_type = param.RK_LSEX; % LS explicit
    LSEXa = r4ka;
    LSEXb = r4kb;
    LSEXc = r4kc;
    num_stages = 5;
end

% flag to tell biot_rhs_3d whether or not to include LF dissipation terms

param.include_lf = ...
    (RK_type == param.RK_LSEX && param.splitting == param.SPL_NONE);

% Storage for residuals, depending on RK scheme. This assumes that
% field and field_adj are correctly initialised; if either is empty,
% so are the corresponding residuals.

if RK_type == param.RK_LSEX
    res_field = zeros(size(field));
    res_mem_hf = zeros(size(mem_hf));
    res_field_adj = zeros(size(field_adj));
    res_mem_hf_adj = zeros(size(mem_hf_adj));
    res_sens_ker = zeros(size(sens_ker));
elseif RK_type == param.RK_LSIMEX
    res_field_IM = zeros(size(field));
    res_field_EX = zeros(size(field));
    res_mem_hf_EX = zeros(size(mem_hf));
    % IMEX for adjoint field not yet implemented
    % IMEX for sensitivity kernels not yet implemented 
end    

% Storage for receiver values. numel(param.elm_rec) is the number of receivers.
% Need num_steps+1 time slots because we store initial values in the first
% place. If we're not simulating either type of receiver, initialise to
% empty (these are return values of biot_3d so must be set to something)

if param.run_forward
    receiver =     zeros(numel(param.elm_rec), num_steps+1, param.Nfields);
else
    receiver = [];
end
if param.run_forward_adj
    receiver_adj = zeros(numel(param.elm_rec), num_steps+1, param.Nfields);
else
    receiver_adj = [];    
end

% Storage for times at RK steps (row vector)

step_time = zeros(1, num_steps+1);

% Initialise sensitivity kernels

if param.run_adjoint_method
    field0 = field; % this is used to compute kernels to estimate time derivative of fields
    ker_rho_a = zeros(size(field,1), size(field,2));
    ker_kappa_fr = zeros(size(field,1), size(field,2));
end


% outer time step loop. Looks like it has one step too many but this is
% just to get the status printout at the end of the loop; the computational
% body of the loop isn't executed when tstep == num_steps

for tstep = 0:num_steps

    time = initial_time + tstep*dt;
    step_time(tstep+1) = time;
    
    % Interpolate field values to give receiver values at this time step
    % for adjoint and non-adjoint wavefields, as appropriate
    
    for j=1:numel(param.elm_rec)
        if param.run_forward
            receiver(j, tstep+1, :) = param.interp_rec(j,:)*squeeze(field(:,param.elm_rec(j),:));
        end
        if param.run_forward_adj
            receiver_adj(j, tstep+1, :) = param.interp_rec(j,:)*squeeze(field_adj(:,param.elm_rec(j),:));
        end        
    end

    % Diagnostics: display some field and timing data every 
    % param.steps_per_report steps and after the final step is completed 
    % (after which the loop is terminated by break statement)
    
    if mod(tstep, param.steps_per_report) == 0 || tstep == num_steps

        fprintf('*************  After %d time step(s) of %d *************\n', tstep, num_steps);
        
        if param.run_forward
            if isfield(param, 'exact')
                exact = param.exact(time, x, y, z, param);
                display_field_2(field, fname, 'Simulated', exact, fname, 'Analytic')
            else % ~isfield(param, 'exact')
                if ~param.run_forward_adj
                    display_field_1(field, fname, 'Simulated')
                end
            end % isfield(param, 'exact')
        end % param.run_forward
        if param.run_forward_adj
            if isfield(param, 'exact_adj')
                exact_adj = param.exact_adj(time, x, y, z, param);
                display_field_2(field_adj, fname_adj, 'Simulated Adjoint', exact_adj, fname_adj, 'Analytic adjoint')
            else % ~isfield(param, 'exact')
                if ~param.run_forward
                    display_field_1(field_adj, fname_adj, 'Simulated Adjoint')
                end
            end % isfield(param, 'exact')
        end % param.run_forward_adj
        if param.run_forward && param.run_forward_adj
            display_field_2(field, fname, 'Simulated', field_adj, fname_adj, 'Simulated Adjoint')
        end
        
        fprintf('Time spent on time integration: %2.3fs [%s] CPU %2.3fs [%s] wallclock\n', cputime-t0, ftime(round(cputime-t0)), toc, ftime(round(toc)));
        if tstep ~= 0
            fprintf('Estimated remaining wallclock time %2.3fs [%s]\n', (num_steps-tstep)/tstep*toc, ftime(round((num_steps-tstep)/tstep*toc)));
        end
        if tstep ==  num_steps
            break
        end
        
    end % mod(tstep, param.steps_per_report) == 0 || tstep ==  num_steps

    % if we're doing operator-splitting with an analytical term first, 
    % advance the low-frequency elements by the analytic solution of the 
    % stiff part of the equation (for half a step, if it's Strang)
    % Operator splitting not yet implemented for adjoint
        
    if param.splitting == param.SPL_SIMPLE_AN 
        field(:,param.lf_elt,:) = analytic_stiff_term(field(:,param.lf_elt,:), param, dt);
    elseif param.splitting == param.SPL_STRANG_ANA
        field(:,param.lf_elt,:) = analytic_stiff_term(field(:,param.lf_elt,:), param, 0.5*dt);
    end
    
    % Select explicit or IMEX Runge-Kutta scheme
    
    if RK_type == param.RK_LSEX
        
        % Initialise to empty so RK updates work even if RHS function
        % not called
        
        rhs_field = [];
        rhs_mem_hf = [];
        rhs_field_adj = []; 
        rhs_mem_hf_adj = [];
        rhs_ker = [];
        
        % This is a low-storage Runge-Kutta method following Carpenter 
        % and Kennedy

        for INTRK = 1:num_stages %  RK stage loop for LSEX
            % time value (used with to compute analytic plane wave solution 
            % and update boundary fluxes)
            RKtime = time+dt*LSEXc(INTRK);
            % compute right hand side. Boolean parameter after time is true
            % for adjoint wavefield calculation, false for non-adjoint
            % wavefield
            if param.run_forward
                [rhs_field, rhs_mem_hf] = biot_rhs_3d(RKtime, false, field, mem_hf, param);
            end
            if param.run_forward_adj % rhs_mem_hf_adj not used at this point
                [rhs_field_adj, rhs_mem_hf_adj] = biot_rhs_3d(RKtime, true, field_adj, mem_hf_adj, param);
            end
            if param.run_adjoint_method
                rhs_ker = kernel_rhs(field, rhs_field, field_adj, param);
            end
            % initiate and increment Runge-Kutta residuals
            res_field     = LSEXa(INTRK)*res_field +     dt*rhs_field;
            res_field_adj = LSEXa(INTRK)*res_field_adj + dt*rhs_field_adj;
            % update fields
            field     = field     + LSEXb(INTRK)*res_field;
            field_adj = field_adj + LSEXb(INTRK)*res_field_adj;
            % attenuation terms (high-frequency case, might be empty)
            % HF dissipation not yet implemented for adjoint
            res_mem_hf = LSEXa(INTRK)*res_mem_hf + dt*rhs_mem_hf;
            mem_hf     = mem_hf + LSEXb(INTRK)*res_mem_hf;
            % sensitivity kernels, might be empty
            res_sens_ker = LSEXa(INTRK)*res_sens_ker + dt*rhs_ker;
            sens_ker     = sens_ker + LSEXb(INTRK)*res_sens_ker;
        end % RK stage loop for LSEX
    elseif RK_type == param.RK_LSIMEX
        % IMEX not yet implemented for adjoint
        % This is the low-storage IMEX method from Cavagliere and Bewley
        % The memory variables have no stiff part so a simplified form
        % is used, in which the implicit register is absent because it
        % always holds 0
        for INTRK = 1:num_stages % RK stage loop for LSIMEX
            RKtime = time+dt*LSIMEXc(INTRK);
            if INTRK == 1
                res_field_EX = field;
                res_mem_hf_EX = mem_hf;
            else % INTRK != 1
                res_field_EX  = field  + (LSIMEXaIMsd(INTRK-1)-LSIMEXb(INTRK-1))*dt*res_field_IM+(LSIMEXaEXsd(INTRK-1)-LSIMEXb(INTRK-1))*dt*res_field_EX;
                res_mem_hf_EX = mem_hf + (LSIMEXaEXsd(INTRK-1)-LSIMEXb(INTRK-1))*dt*res_mem_hf_EX;
            end % if INTRK == 1
            res_field_IM(:,param.lf_elt,:) = IMEX_stiff_term(LSIMEXaIMd(INTRK)*dt, res_field_EX(:,param.lf_elt,:), param);
            [res_field_EX, res_mem_hf_EX] = biot_rhs_3d(...
                RKtime, ...
                res_field_EX + LSIMEXaIMd(INTRK)*dt*res_field_IM, ...
                res_mem_hf_EX, param);
            field  = field  + LSIMEXb(INTRK)*dt*res_field_IM+LSIMEXb(INTRK)*dt*res_field_EX;
            mem_hf = mem_hf + LSIMEXb(INTRK)*dt*res_mem_hf_EX;
        end % RK stage loop for LSIMEX
    else % unknown Rk type
        error('RK_type %d not recognised', RK_type)
    end % if RK_type == 
    
    % if we're doing operator-splitting with an analytical term second, 
    % advance the low-frequency elements by the analytic solution of the 
    % stiff part of the equation (for half a step, if it's Strang)
    % Operator splitting not yet implemented for adjoint
    
    if param.splitting == param.SPL_SIMPLE_NA 
        field(:,param.lf_elt,:) = analytic_stiff_term(field(:,param.lf_elt,:), param, dt);
    elseif param.splitting == param.SPL_STRANG_ANA
        field(:,param.lf_elt,:) = analytic_stiff_term(field(:,param.lf_elt,:), param, 0.5*dt);
    end

    % Compute sensibility kernels. TODO Need to implement other kernels as per
    % paper
    if param.run_adjoint_method
        [ker_rho_a, ker_kappa_fr] = compute_kernels(ker_rho_a, ker_kappa_fr, dt, field, field0, field_adj, param);
        % store previous time step temporarily for computation of ker_rho_a
        field0 = field;
    end
    
end % time-stepping loop

% correct last time step to initial_time + time_width in case 
% initial_time + num_steps*dt differs slightly because of FP noise

step_time(end) = initial_time + time_width;

% Package all output data into single struct to return

output.step_time = step_time;
output.field = field;
output.mem_hf = mem_hf;
output.receiver = receiver;
output.field_adj = field_adj;
output.mem_hf_adj = mem_hf_adj;
output.receiver_adj = receiver_adj;
output.ker_rho_a = ker_rho_a;
output.ker_kappa_fr = ker_kappa_fr;
output.sens_ker = sens_ker;

return

% Utilities for displaying field magnitudes and differences

function display_field_1(fieldX, nameX, titleX)

fprintf('%s\n', titleX);
fprintf('Field Max\n');
for fld = 1:numel(nameX)
    simulX = fieldX(:, :, fld);
    simulX = simulX(:);
    fprintf('%-5s %.7e\n', nameX{fld}, max(abs(simulX)));
end

function display_field_2(fieldX, nameX, titleX, fieldY, nameY, titleY)

fprintf('%-20s%-20s\n', titleX, titleY);
fprintf('Field Max           Field Max           Max diff abs  Max diff rel\n')
for fld = 1:numel(nameX)
    simulX = fieldX(:, :, fld);
    simulX = simulX(:);
    mX = max(abs(simulX));
    simulY = fieldY(:, :, fld);
    simulY = simulY(:);
    mY = max(abs(simulY));
    mXY = max(abs(simulX-simulY));
    fprintf('%-5s %.7e %-5s %.7e', nameX{fld}, mX, nameY{fld}, mY);
    if strcmp(nameX{fld}, nameY{fld})
        fprintf(' %.7e', mXY);
        if max(abs(simulY)) > 0
            fprintf(' %.7e', mXY/mY);
        end
    end
    fprintf('\n');
end
