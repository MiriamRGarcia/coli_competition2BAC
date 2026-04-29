%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CostFun_Fit_FCdata: Cost function (unweighted least squares) to fit the
% model against flow cytometry data for the competition experiment 
% between strains:
% Wild-type E. coli YFP (W-YFP) vs Tolerant E. coli mCherry (T-mCh),
% Wild-type E. coli mCherry (W-mCh) vs Tolerant E. coli YFP (T-YFP).
% The competing strains were exposed to sucessive cycles alternating 
% growth periods and disinfection periods with benzalkonium chloride (BAC). 
% After each disinfection step, the cultures were diluted into fresh M9
% medium to start the next treatment cycle (4 complete cycles).
% Growth periods lasted 24 hours.
% Disinfection periods lasted 10 minutes.
% Flow cytometry (FC) data was collected before each disinfection round.
% Cell count data were collected before each disinfection round.
% The model was fitted to FC data to compare the fitted
% microbial traits with those determined from independent experiments
% without differentiating cells labelled with YFP and mCherry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = CostFun_Fit_FCdata(decVars, input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% (1) Initialise variables:

nc     = input.nc;
BAC    = input.BAC;
X_W0   = input.X_W0;
X_T0   = input.X_T0;
X_e    = input.X_e;
Data_W = input.Data_W;
Data_T = input.Data_T;
    
% Problem sizes:
nBAC   = numel(BAC);
nRe    = size(Data_W, 4);
nT0_W0 = numel(X_W0);

% Obtain parameter values:
mu_W = decVars(1);
mu_T = decVars(nBAC + 1);

if size(decVars, 1) > 1
    SF_W  = [1;decVars(2:nBAC)];
    SF_T  = [1;decVars((nBAC+2):end)];
else
    SF_W  = [1;decVars(2:nBAC).'];
    SF_T  = [1;decVars((nBAC+2):end).'];
end

% ----------------------------------------------------------------------- %
% (2) Simulate model:

Model_W = zeros(nc, nBAC, nT0_W0);
Model_T = Model_W;

for iT0_W0 = 1:nT0_W0

    % Initial condition:
    x_0 = [X_W0(iT0_W0);X_T0(iT0_W0)];

    % Loop in the BAC concentrations:
    for iBAC = 1:nBAC
        
        % Call the simulation function:
        input.mu_W = mu_W;
        input.mu_T = mu_T;
        input.SF_W = SF_W(iBAC);
        input.SF_T = SF_T(iBAC);
        input.x_0  = x_0;

        output     = Sim_Model(input);
        X_W_bef    = output.X_W_bef;
        X_T_bef    = output.X_T_bef;

        % Save model:
        ext_ind                        = find(max(X_W_bef, X_T_bef) < X_e);
        X                              = X_W_bef + X_T_bef;
        Model_W(1:nc, iBAC, iT0_W0)    = X_W_bef./X; 
        Model_T(1:nc, iBAC, iT0_W0)    = X_T_bef./X;
        Model_W(ext_ind, iBAC, iT0_W0) = 0;
        Model_T(ext_ind, iBAC, iT0_W0) = 0;
    end
end

% ----------------------------------------------------------------------- %
% (3) Calculate cost function:

Model_W = reshape(repmat(Model_W, 1, 1, 1, nRe), [], 1);
Model_T = reshape(repmat(Model_T, 1, 1, 1, nRe), [], 1);

Data_W  = reshape(Data_W, [], 1);
Data_T  = reshape(Data_T, [], 1);

J = sum((Model_W - Data_W).^2) + sum((Model_T - Data_T).^2);

end