%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sim_Model: Function to simulate model of the competition experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = Sim_Model(input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% (1) Initialise variables:

% Obtain setup parameters:
nc   = input.nc;
t_g  = input.t_g;
tsim = input.tsim;
D    = input.D;
K    = input.K;
X_e  = input.X_e;
x_0  = input.x_0;

% Obtain microbial traits:
mu_W = input.mu_W;                                   
mu_T = input.mu_T;                                                                                              
SF_W = input.SF_W;                                                     
SF_T = input.SF_T;                                                           

% Initialise saturation times:
tsat = zeros(nc, 1);

% Initialise cell densities:
X_W         = [];
X_T         = [];

% Initialise cell densities before killing:
X_W_bef     = [];
X_T_bef     = [];

% Initialise flags to detect population extinction:
if x_0(1) < X_e,ext_flag_W = 1;else,ext_flag_W = 0;end
if x_0(2) < X_e,ext_flag_T = 1;else,ext_flag_T = 0;end

% ----------------------------------------------------------------------- %
% (2) Simulate selection experiment:

for ic = 1:nc
    
    % Obtain simulation times within the current growth period:
    tsim_aux = tsim(tsim >= (ic - 1)*t_g & tsim <= ic*t_g);

    % Initial time of the current period:
    t0 = tsim_aux(1);
    
    % Simulate exponential growth within the current period:
    if ext_flag_W > 0
        X_W_ic = x_0(1)*ones(size(tsim_aux));
    else
        X_W_ic = x_0(1)*exp(mu_W*(tsim_aux - t0));
    end
    if ext_flag_T > 0
        X_T_ic = x_0(2)*ones(size(tsim_aux));
    else
        X_T_ic = x_0(2)*exp(mu_T*(tsim_aux - t0));
    end
    
    % Total population in the growth period:
    X_ic  = X_W_ic + X_T_ic;
    
    % Find saturation time of the growth period:
    ind_tsat = find(X_ic <= K, 1, 'last');
    tsat(ic) = tsim_aux(ind_tsat) - t0;
    
    % Rewrite states during the growth period for the stationary phase:
    % (note! numel(indtR)>= 1 when the total initial concentration is
    % <= K, which is always fullfiled. If the stationary phase is not
    % reached, the command below has no effect):
    X_W_ic(ind_tsat:end) = X_W_ic(ind_tsat);
    X_T_ic(ind_tsat:end) = X_T_ic(ind_tsat);

    % Actualise states:
    X_W = [X_W;X_W_ic(1:end-1)];
    X_T = [X_T;X_T_ic(1:end-1)];
    
    % Actualise cell densities before killing:
    X_W_bef = [X_W_bef;X_W_ic(end)];
    X_T_bef = [X_T_bef;X_T_ic(end)];
    
    % Actualise initial condition (instantaneous killing):
    X_W0 = D*SF_W*X_W_ic(end);
    X_T0 = D*SF_T*X_T_ic(end);  
        
    % Extinction conditions:       
    if X_W0 < X_e, ext_flag_W = 1;end
    if X_T0 < X_e, ext_flag_T = 1;end
    
    % Actualise initial condition:
    x_0 = [X_W0;X_T0];

    % Save the last simulation time:
    if ic > nc - 1
        X_W = [X_W;X_W0];
        X_T = [X_T;X_T0];
   end
end

% Obtain output variables:
output.X_W     = X_W;
output.X_T     = X_T;
output.X_W_bef = X_W_bef;
output.X_T_bef = X_T_bef;
output.tsat    = tsat;

end
