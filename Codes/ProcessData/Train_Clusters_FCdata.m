%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train_Clusters_FCdata: Train Gaussian Mixture Model (GMM) clustering 
% with a subset of wells from flow cytometry (FC) data of the 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables

% Set seed to generate randoms for reproducibility:
rng('default')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% Wells to calibrate YFP/mCherry signals from W and T:
wells_YFP = {'E4', 'E5', 'E6', 'F4', 'F5', 'F6'}; % Wells with mostly YFP (day 1);
                         
wells_mCh = {'E7', 'E8', 'E9', 'F7', 'F8', 'F9'}; % Wells with only mCherry (day 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Load data of the wells for training the clusters:

% Number of wells used for calibration of YFP:
nw     = size(wells_YFP, 2);
YFP_FL = [];

for iw = 1:nw
    
    % Read .mat file:
    load('Data_CompExp_FC.mat', wells_YFP{iw})

    % Signal of YFP on FITC and PE channels:
    eval(sprintf('aux = [%s.FL1_H{1} %s.FL5_H{1}];', wells_YFP{iw}, wells_YFP{iw}))
    YFP_FL = [YFP_FL;aux];   
end

% Number of wells used for calibration of mCherry:
nw     = size(wells_mCh, 2);
mCh_FL = [];

for iw = 1:nw
    
    % Read .mat file:
    load('Data_CompExp_FC.mat', wells_mCh{iw})

    % Signal of mCherry on FITC and PE channels:
    eval(sprintf('aux = [%s.FL1_H{1} %s.FL5_H{1}];', wells_mCh{iw}, wells_mCh{iw}))
    mCh_FL = [mCh_FL;aux];
end

% ----------------------------------------------------------------------- %
%% (2) Fit the mixed Gaussian to the arquetype signals:

% Options for fitgmdist:
opts    = statset('MaxIter', 1500, 'TolFun', 1e-12);

% Fit the Gaussian mixture model model to data:
Z   = [log10(1 + YFP_FL);log10(1 + mCh_FL)];
GMM = fitgmdist(Z , 3, 'Options', opts);


% ----------------------------------------------------------------------- %
%% (3) Classify data in clusters:

% Size of the dataset:
nData   = size(Z, 1);

% Posterior probabilities of component-member membership:
P       = posterior(GMM, Z);

% Set the interval of probabilities to be member of both clusters:
P_thold = [0.45 0.55];

% Clasification of the points in clusters according to the mixture model:
idClus  = cluster(GMM, Z);

% Indexes of the points belonging to both clusters:
idBothClus12  = find(P(:,1) >= P_thold(1) & P(:,1) <= P_thold(2)); 
idBothClus13  = find(P(:,2) >= P_thold(1) & P(:,2) <= P_thold(2)); 
idBothClus23  = find(P(:,3) >= P_thold(1) & P(:,3) <= P_thold(2)); 

% Obtain indexes of the points in the sorted array:
[~, ord] = sort(P(:,1));

% ----------------------------------------------------------------------- %
%% (4) Plot clustering results:

% Colours for clusters:
cols = [229 127 127;
        38 91 143;
        207 179 65]/256;
    
% Plot cluster classification and contours:
figure
hold on
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
gscatter(Z(:, 1), Z(:, 2), idClus, cols, '+ox', 5)
plot(Z(idBothClus12,1), Z(idBothClus12,2), 'ko', 'MarkerSize', 10)
plot(Z(idBothClus23,1), Z(idBothClus23,2), 'ko', 'MarkerSize', 10)

% Contours:
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMM,[x0,y0]),x,y);
fcontour(gmPDF, [2.5 4 2 4])
title('Fitted GMM clusters and contours', 'Interpreter', 'Latex', 'FontSize', 12)
xlabel('FITC-H', 'Interpreter', 'Latex', 'FontSize', 12)
ylabel('PE-H', 'Interpreter', 'Latex', 'FontSize', 12)
hold off
legend off
xlim([2 4])
ylim([1 4.5])

% ----------------------------------------------------------------------- %
% Save results:
save('Data_CompExp_TrainedGMM', 'GMM', 'P_thold')