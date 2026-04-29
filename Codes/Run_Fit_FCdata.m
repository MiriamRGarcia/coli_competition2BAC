%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_Model: Run model model fitting against flow cytometry data for the
% competition experiment between:
% Wild-type E. coli YFP (W-YFP) vs Tolerant E. coli mCherry (T-mCh),
% Wild-type E. coli mCherry (W-mCh) vs Tolerant E. coli YFP (T-YFP).
% The competing strains were exposed to sucessive cycles alternating 
% growth periods and disinfection periods with benzalkonium chloride (BAC). 
% After each disinfection step, the cultures were diluted into fresh M9
% medium to start the next treatment cycle (4 complete cycles).
% Growth periods lasted 24 hours.
% Disinfection periods lasted 10 minutes.
% Flow cytometry (FC) data was collected before each disinfection round.
% The model was fitted to FC data to compare the fitted
% microbial traits with those determined from independent experiments
% without differentiating cells labelled with YFP and mCherry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% Duration of the growth periods (min):
t_g       = 24*60;

% Time step for model simulation (min):
ht        = 1/60;

% BAC doses of the competition experiment (mug/ml):
BAC       = [0;30;40;50];

% Maximum optimisation time for model fitting (s):
max_topt  = 10*60;

% Carrying capacity (CFUS/mL):
K         = 5e8;

% Dilution factor:
D         = 1/100;

% Total initial inocula:
X_0      = 1e6;

% Initial mixing ratios (X_T0/X_W0):
T0_W0    = [1;1e-2;1e-4];

% Culture volume (mL):
V        = 600e-3;

% Extinction limit (CFU/mL):
X_e      = 1/V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Initialise variables:

% Initialise names of the variables:
strains     = {'W_YFP', 'W_mCh', 'T_YFP', 'T_mCh'};
c_names_W   = {'W_YFP', 'W_mCh'};
c_names_T   = {'T_mCh', 'T_YFP'};
opt_names   = {'opt_pars_W_YFP', 'opt_pars_W_mCh'};

% Number of competition cases:
n_comp = size(c_names_W, 2);

% Confidence bounds for model fit:
load('ProcessData/Data_Traits.mat', 'CI_mu_W', 'CI_mu_T', 'CI_SF_W', 'CI_SF_T')

% Load flow cytometry (FC) data for the different cases:
load('ProcessData/Data_CompExp_FC_Clustered.mat', 'W_YFP_1_0', 'W_YFP_1_2', 'W_YFP_1_4',...
     'W_mCh_1_0', 'W_mCh_1_2', 'W_mCh_1_4');

% Remove bad replicate:
W_mCh_1_0.mCh_freq(1, 1, 1) = W_mCh_1_0.mCh_freq(1, 1, 2);
W_mCh_1_0.YFP_freq(1, 1, 1) = W_mCh_1_0.YFP_freq(1, 1, 2);

% Preprocess data:
nc      = size(W_mCh_1_0.YFP_freq, 1);
nBAC    = size(W_mCh_1_0.YFP_freq, 2);
nRe     = size(W_mCh_1_0.YFP_freq, 3);
nT0_W0  = numel(T0_W0);

Data_W_YFP = zeros(nc, nBAC, nT0_W0, nRe);
Data_W_mCh = zeros(nc, nBAC, nT0_W0, nRe);
Data_T_YFP = zeros(nc, nBAC, nT0_W0, nRe);
Data_T_mCh = zeros(nc, nBAC, nT0_W0, nRe);

Data_W_YFP(1:nc, 1:nBAC, 1, 1:nRe) = W_YFP_1_0.YFP_freq;
Data_W_YFP(1:nc, 1:nBAC, 2, 1:nRe) = W_YFP_1_2.YFP_freq;
Data_W_YFP(1:nc, 1:nBAC, 3, 1:nRe) = W_YFP_1_4.YFP_freq;
Data_W_mCh(1:nc, 1:nBAC, 1, 1:nRe) = W_mCh_1_0.mCh_freq;
Data_W_mCh(1:nc, 1:nBAC, 2, 1:nRe) = W_mCh_1_2.mCh_freq;
Data_W_mCh(1:nc, 1:nBAC, 3, 1:nRe) = W_mCh_1_4.mCh_freq;

Data_T_mCh(1:nc, 1:nBAC, 1, 1:nRe) = W_YFP_1_0.mCh_freq;
Data_T_mCh(1:nc, 1:nBAC, 2, 1:nRe) = W_YFP_1_2.mCh_freq;
Data_T_mCh(1:nc, 1:nBAC, 3, 1:nRe) = W_YFP_1_4.mCh_freq;
Data_T_YFP(1:nc, 1:nBAC, 1, 1:nRe) = W_mCh_1_0.YFP_freq;
Data_T_YFP(1:nc, 1:nBAC, 2, 1:nRe) = W_mCh_1_2.YFP_freq;
Data_T_YFP(1:nc, 1:nBAC, 3, 1:nRe) = W_mCh_1_4.YFP_freq;


% Initial concentrations of W and S4:
X_W0 = X_0./(1 + T0_W0);
X_T0 = X_0 - X_W0; 

% Simulation times (min):
tsim = (0:ht:nc*t_g).';

% ----------------------------------------------------------------------- %
%% (2) Set ESS options:

% General options:
opts.maxtime      = max_topt;                                                                           
opts.maxeval      = 1.0e10;                              
opts.strategy     = 3;                                        
opts.local.solver = 'fminsearch';                                                                            
opts.local.finish = 'fminsearch';                                      
opts.local.n1     = 5;                                                    
opts.local.n2     = 5;

% Name of the cost function:
problem.f = 'CostFun_Fit_FCdata';

%- ---------------------------------------------------------------------- %
%% (3) Fit the model to FC data:

% Set lower bounds:
problem.x_L = [CI_mu_W(1);CI_SF_W(:, 1);CI_mu_T(1);CI_SF_T(:, 1)];

% Set upper bounds:
problem.x_U = [CI_mu_W(2);CI_SF_W(:, 2);CI_mu_T(2);CI_SF_T(:, 2)];

% Set initial condition:
problem.x_0 = problem.x_L + (problem.x_U - problem.x_L).*rand(size(problem.x_L));
    
% Perform parameter estimation for each case:
for i_comp = 1:n_comp

    % Main call to ESS:
    eval(sprintf('auxData_W = Data_%s;', c_names_W{i_comp}))
    eval(sprintf('auxData_T = Data_%s;', c_names_T{i_comp}))

    % Remove NaN values when total population is extinct:
    auxData_W(isnan(auxData_W)) = 0;
    auxData_T(isnan(auxData_T)) = 0;

    % Call to optimisation function (ESS):
    input.nc     = nc;
    input.t_g    = t_g;
    input.tsim   = tsim;
    input.BAC    = BAC;
    input.D      = D;
    input.K      = K;
    input.X_W0   = X_W0;
    input.X_T0   = X_T0;
    input.X_e    = X_e;
    input.Data_W = auxData_W;
    input.Data_T = auxData_T;
    
    Results = ess_kernel(problem, opts, input);

    % Obtain optimal parameters:
    opt_pars = Results.xbest.';
    
    eval(sprintf('%s = opt_pars;', opt_names{i_comp}))

    % Remove mat file generated by ESS:
    delete ess_report.mat
    
    % Print results:
    fprintf('\n >> The fitted parameters for the competition case: %s vs %s are:\n', c_names_W{i_comp}, c_names_T{i_comp})
    fprintf('>> Growth rate of strain %s (1/h): %.2f\n', c_names_W{i_comp}, 60*opt_pars(1))
    fprintf('>> Growth rate of strain %s (1/h): %.2f\n', c_names_T{i_comp}, 60*opt_pars(nBAC+1))
    for iBAC = 2:nBAC
        fprintf('>> Survival fraction of strain %s at BAC = %u: %.02e\n', c_names_W{i_comp}, BAC(iBAC), opt_pars(iBAC))
        fprintf('>> Survival fraction of strain %s at BAC = %u: %.02e\n', c_names_T{i_comp}, BAC(iBAC), opt_pars(nBAC+iBAC))
    end
end


%- ---------------------------------------------------------------------- %
%% (4) Plot results:

% Define colors for model predictions:
col_W    = [206 92 92]/256;
col_T    = [183 145 47]/256;

% Define colors for data replicates:
col_Re_W = [255 209 209;
            255 182 182;
            255 142 142;
            229 127 127;
            206 114 114;
            185 102 102]/256;
        
col_Re_T = [241 223 152;
            235 210 109;
            207 179 65;
            186 161 58;
            167 144 52;
            133 115 41]/256;
        
% Data sizes:
nc_data    = nc;

% Values of the x-axis to plot model:
hplot      = 2;
tplot      = (1:hplot:nc*hplot).';

% Labels for the x-axis (cycles):
tlabs      = (1:nc).';

% Values of the x-axis to plot data:
tplot_data = (1:hplot:nc_data*hplot).';

input.tplot      = tplot;
input.tplot_data = tplot_data;
input.tlabs      = tlabs;
input.BAC        = BAC;
input.col_W      = col_W;
input.col_T      = col_T;
input.col_Re_W   = col_Re_W;
input.col_Re_T   = col_Re_T;

for i_comp = 1:n_comp
    
    % Names of the competition cases:
    if n_comp > 1
        fig_name = strcat('W-',c_names_W{i_comp}(4:end),'\&','T-', c_names_T{i_comp}(4:end));
    else
        fig_name = 'W \& T';
    end

    % Initialise model:
    X_W_bef = zeros(nc, nBAC, nT0_W0);
    X_T_bef = X_W_bef;

    % Simulate model with the optimal parameters:
    eval(sprintf('opt_pars = %s;', opt_names{i_comp}))

    % Obtain parameter values:
    mu_W = opt_pars(1);
    mu_T = opt_pars(nBAC+1);
    SF_W = [1;opt_pars(2:nBAC)];
    SF_T = [1;opt_pars((nBAC+2):end)];

    for iT0_W0 = 1:nT0_W0

        for iBAC = 1:nBAC
            
            % Call the simulation function:
            input.mu_W = mu_W;
            input.mu_T = mu_T;
            input.SF_W = SF_W(iBAC);
            input.SF_T = SF_T(iBAC);
            input.x_0  = [X_W0(iT0_W0);X_T0(iT0_W0)];

            output     = Sim_Model(input);
            
            X_W_bef(1:nc, iBAC, iT0_W0) = output.X_W_bef;
            X_T_bef(1:nc, iBAC, iT0_W0) = output.X_T_bef;
        end
    end

    % Call plot function:
    eval(sprintf('auxData_W = Data_%s;', c_names_W{i_comp}))
    eval(sprintf('auxData_T = Data_%s;', c_names_T{i_comp}))

    % Remove NaN values when total population is extinct:
    auxData_W(isnan(auxData_W)) = 0;
    auxData_T(isnan(auxData_T)) = 0;
    
    % Call plot function:
    input.fig_name  = fig_name;
    input.Model_W   = X_W_bef;
    input.Model_T   = X_T_bef;
    input.FC_Data_W = auxData_W;
    input.FC_Data_T = auxData_T;
    Plot_Model_vs_Data(input)
end


% Save results:
save('ProcessData/Data_CompExp_Fitted_Traits.mat')

rmpath('Functions')
