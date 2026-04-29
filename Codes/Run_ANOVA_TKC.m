%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_ANOVA_TKC: Run two-way ANOVA with repeated measures on data from 
% the preliminary time-kill assay to discern if fluorescent labelling 
% (YFP or mCherry) affects cell survival.
% Strains:
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

% ----------------------------------------------------------------------- %
%% (1) Initialise variables:

% Strains:
strains = {'W_YFP', 'W_mCh', 'T_YFP', 'T_mCh'};

% Load data from the time-kill assay:
load('ProcessData/Data_TKC.mat', 'nRe', 'BAC', 'ts', strains{:})

% BAC concentrations for the competition experiment:
BAC_CE = [30 40 50].';
indBAC = find(ismember(BAC, BAC_CE));
BAC    = BAC_CE;

% Problem sizes:
nBAC = numel(BAC);
nts  = numel(ts);
nstr = size(strains, 2);

% ----------------------------------------------------------------------- %
%% (2) Perform ANOVA globally (with within-subject variables):

% ----------------------------------------------------------------------- %
% (2.a) Construct table with data:

% Define values for labels:
StrainLabel = [repmat('YFP', nRe, 1);repmat('mCh', nRe, 1);repmat('YFP', nRe, 1);repmat('mCh', nRe, 1)];

% Define values for strain name:
StrainName  = [repmat('W', nRe, 1);repmat('W', nRe, 1);repmat('T', nRe, 1);repmat('T', nRe, 1)];

% Define names for the head of the datatable:
ColsNames    = cell(1, 2 + nBAC*nts);
ColsNames{1} = 'Strain';
ColsNames{2} = 'Label';
for iBAC = 1:nBAC
    for its = 1:nts
        ColsNames{2 + (iBAC - 1)*nts + its} = sprintf('BAC%s_t%s', num2str(BAC(iBAC)), num2str(ts(its)));
    end
end

% Resize TKC to an appropiate format:
data = [];
for istr = 1:nstr
    
    % Obtain data in log-scale:
    eval(sprintf('data_aux =  log10(max(1, %s(1:nts, 1:nBAC, 1:nRe)));', strains{istr}))
    
    % Reshape data array to have replicates in rows, times in colums and
    % BAC concentrations in the third dimension:
    data_aux = permute(data_aux, [3 1 2]);
    
    % Reshape data array to have replicates in rows, and the different
    % combinations of times and BAC in columns (ordering = BAC30-t0,
    % BAC-30-t5, BAC30-t10, BAC30-t20, ..., BAC60-t0, BAC60-t5, BAC60-t10,
    % BAC60-t20):
    data_aux = reshape(data_aux, nRe, []);
    
    % Concatenate data for the different cases:
    data = [data;data_aux];
end

% Construct datatable:
between = table(StrainName, StrainLabel, data(:,1), data(:,2), data(:,3), data(:,4),...
          data(:, 5), data(:, 6), data(:, 7), data(:, 8), data(:, 9),...
          data(:, 10), data(:, 11), data(:, 12), 'VariableNames', ColsNames);

% ----------------------------------------------------------------------- %      
% (2.b) Fit regression model for ANOVA:

% Define within subject variables (times and BAC): 
within = table(repmat(ts, nBAC, 1), reshape(repmat(BAC.', nts, 1), [], 1),...
               'VariableNames', {'Time', 'BAC'}); 
           
rm = fitrm(between, sprintf('%s-%s ~ Strain*Label', ColsNames{3}, ColsNames{end}), 'WithinDesign', within);

% ----------------------------------------------------------------------- %
% (2.c) Perform ANOVA:
ranovatbl = ranova(rm, 'WithinModel', 'Time+BAC')
