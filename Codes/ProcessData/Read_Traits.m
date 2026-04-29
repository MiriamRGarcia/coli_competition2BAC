%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read_Traits: Obtain experimental values and confidence intervals (CIs)
% for the microbial traits (growth rates and log(survival fraction)) for 
% the strains of the competition experiment,assuming Gaussian sample of
% the microbial traits.
% The competing strains were exposed to sucessive cycles alternating 
% growth periods and disinfection periods with benzalkonium chloride. 
% After each disinfection step, the cultures were diluted into fresh M9
% medium to start the next cycle.
% Growth periods lasted 24 hours.
% Disinfection periods lasted 10 minutes.
% Growth rates for the competition experiment were obtained from a previous
% study for strains W and T without fluorescent labels.
% Survival fractions were obtained from time-kill curves (TKC) after
% 10 minutes of disinfection for strains:
% Wild-type E. coli tagged with YFP (W-YFP),
% Wild-type E. coli tagged with mCherry (W-mCh),
% Tolerant E. coli tagged with YFP (T-YFP),
% Tolerant E. coli tagged with mCherry (T-mCh).
% TKCs were obtained from an overnight culture with inocula 1e5 CFU/ml,
% following growth during 24h in M9 medium to reach the stationary phase, 
% and then applying benzalkonium chloride (BAC) at concentrations:
% BAC = 30,40,50,60,75,150 microg/ml. 
% The sampling times of the TKC were t = 0,5,10,20 min.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

% Set seed to generate randoms for reproducibility:
rng('default')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% (nonzero) BAC concentrations (microg/mL) for the competition experiment:
BAC_CE    = [30 40 50].';

% Low detection limit (CFU/mL):
LDL       = 100;

% Confidence level for CIs:
confLev   = 0.95;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Initialise variables:

% Confidence variables:
alph  = 1 - confLev;
pLo   = alph/2;
pUp   = 1 - alph/2;

% ----------------------------------------------------------------------- %
%% (2) Confidence intervals for the growth rates:

% Names for the growth rates:
strains  = {'W', 'T'};

% Number of strains:
nstr     = size(strains, 2);

% Sample of the growth rates (1/min):
mu_W     = [1.412597149;1.084390141;1.247790188;1.357616355;1.108894778;1.238256585;1.151457126;1.183867246;1.169272782;1.075031235]/60;
mu_T     = [1.114060307;0.885614791;1.081966982;0.98003761;0.841118003;0.984211572;0.987961864;1.053057416;0.884861726;0.846158344]/60;

% Average growth rates:
ave_mu_W = mean(mu_W);
ave_mu_T = mean(mu_T);

for istr = 1:nstr

    % Obtain sample of growth rates for the strain:
    eval(sprintf('sample = mu_%s;',  strains{istr}))
    
    % Size of the sample:
    n_sample = numel(sample);

    % Estimator of the average:
    eval(sprintf('thetahat = ave_mu_%s;', strains{istr}))
    
    % Standard deviation:
    thetastd = std(sample); 
    
    % Calculate confidence interval for the mean cell counts:
    ndf      = n_sample - 1;
    tStud    = tinv([pLo pUp], ndf);
    CI       = (thetahat + tStud*thetastd/sqrt(n_sample)).';
    
    % CI for the growth rate:
    eval(sprintf('CI_mu_%s = CI;', strains{istr}))
end

% Save average growth rates and confidence intervals:
save('Data_Traits.mat', 'ave_mu_W', 'ave_mu_T', 'mu_W', 'mu_T', 'CI_mu_W', 'CI_mu_T')

% ----------------------------------------------------------------------- %
%% (3) Confidence intervals for the survival fractions (log scale):

% ----------------------------------------------------------------------- %
% (3.a) Differentiating labels:

% Names of the different cases:
strains      = {'W_YFP', 'W_mCh', 'T_YFP', 'T_mCh'};
 
% Names of the survival fractions:
SF_names     = {'SF_W_YFP', 'SF_W_mCh', 'SF_T_YFP', 'SF_T_mCh'};

% Names of the CIs for the survival fractions:
ave_SF_names = {'ave_SF_W_YFP', 'ave_SF_W_mCh', 'ave_SF_T_YFP', 'ave_SF_T_mCh'};

% Names of the CIs for the survival fractions:
CI_SF_names  = {'CI_SF_W_YFP', 'CI_SF_W_mCh', 'CI_SF_T_YFP', 'CI_SF_T_mCh'};

% Number of cases:
nstr         = size(strains, 2);

% Load data from time-kill assay:
load('Data_TKC.mat', 'BAC', SF_names{:})

% Find BAC concentrations for the competition experiment within those
% of the time-kill assay:
indBAC = find(ismember(BAC, BAC_CE));
BAC    = BAC_CE;
nBAC   = numel(BAC);

for istr = 1:nstr
    
    % Obtain replicates of the survival fraction:
    eval(sprintf('SF = %s(indBAC, :);', SF_names{istr}))  
        
    for iBAC = 1:nBAC

        % Obtain replicates of the survival fraction (log) for the BAC conc:
        sample = log(reshape(SF(iBAC, :), [], 1));

        % Remove NaN values (bad replicate for W-mCh):
        NaN_ind       = isnan(sample);
        sample(NaN_ind) = [];

        % Size of the sample:
        n_sample = numel(sample);

        % Estimator of the average normal:
        thetahat = mean(sample);
        
        % Standard deviation:
        thetastd = std(sample);

        % Calculate confidence interval for the mean:
        ndf    = n_sample - 1;
        tStud  = tinv([pLo pUp], ndf);
        CI     = (thetahat + tStud*thetastd/sqrt(n_sample)).';

        % Save  average of the survival fraction:
        eval(sprintf('%s(iBAC, 1) = exp(thetahat);', ave_SF_names{istr}))
        
        % Save  CI for the survival fraction:
        eval(sprintf('%s(iBAC, 1:2) = exp(CI);', CI_SF_names{istr}))
    end
end

% Save average growth rates and confidence intervals:
save('Data_Traits.mat', ave_SF_names{:}, CI_SF_names{:}, "-append")

% ----------------------------------------------------------------------- %
% (3.b) Without differentiating labels:

% Names of the different cases:
strains      = {'W', 'T'};
 
% Names of the survival fractions:
SF_names     = {'SF_W', 'SF_T'};

% Names of the CIs for the survival fractions:
ave_SF_names = {'ave_SF_W', 'ave_SF_T'};

% Names of the CIs for the survival fractions:
CI_SF_names  = {'CI_SF_W', 'CI_SF_T'};

% Number of cases:
nstr         = size(strains, 2);

% Load data from time-kill assay:
load('Data_TKC.mat', 'BAC', SF_names{:})

for istr = 1:nstr
    
    % Obtain replicates of the survival fraction:
    eval(sprintf('SF = %s(indBAC, :);', SF_names{istr}))  
        
    for iBAC = 1:nBAC

        % Obtain replicates of the survival fraction (log) for the BAC conc:
        sample = log(reshape(SF(iBAC, :), [], 1));

        % Remove NaN values (bad replicate for W-mCh):
        NaN_ind       = isnan(sample);
        sample(NaN_ind) = [];

        % Size of the sample:
        n_sample = numel(sample);

        % Estimator of the average normal:
        thetahat = mean(sample);
        
        % Standard deviation:
        thetastd = std(sample);

        % Calculate confidence interval for the mean:
        ndf    = n_sample - 1;
        tStud  = tinv([pLo pUp], ndf);
        CI     = (thetahat + tStud*thetastd/sqrt(n_sample)).';

        % Save  average of the survival fraction:
        eval(sprintf('%s(iBAC, 1) = exp(thetahat);', ave_SF_names{istr}))
        
        % Save  CI for the survival fraction:
        eval(sprintf('%s(iBAC, 1:2) = exp(CI);', CI_SF_names{istr}))
    end
end

% Save average growth rates and confidence intervals:
save('Data_Traits.mat', ave_SF_names{:}, CI_SF_names{:}, "-append")
