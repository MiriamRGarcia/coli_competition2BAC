%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read_TKC: Read data from time-kill curves (TKC) of the preliminary 
% time-kill assay performed with strains:
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% Number of replicates:
nRe       = 4;

% Sampling times (min):
ts        = [0 5 10 20].';

% BAC concentrations (mug/ml):
BAC       = [30 40 50 60 75 150].';

% Low detection limit (CFU/mL):
LDL       = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Read data from Excel file:
filename  = '../../Data/PreliminaryTKC_WvsT_YFP_mCherry_BAC_M9.xlsx';
sheetname = 'TKC';

% Initialise row counter:
irow_W_YFP = 2;
irow_W_mCh = 26;
irow_T_YFP = 50;
irow_T_mCh = 74;

% Rows for the first data of each time:
irow_ts_W_YFP = [98;194;290;0];
irow_ts_W_mCh = [122;218;314;0];
irow_ts_M_YFP = [146;242;338;0];
irow_ts_M_mCh = [170;266;362;0];

% Initialise cell counts, where:
% 1st element: Time;
% 2nd element: BAC concentration;
% 3th element: Replicate;
nts  = numel(ts);
nBAC = numel(BAC);

W_YFP = zeros(nts, nBAC, nRe);
W_mCh = W_YFP;
T_YFP = W_YFP;
T_mCh = W_YFP;

% Loop in the sampling times:
for its = 1:nts
    
    % Loop in the BAC concentrations:
    for iBAC = 1:nBAC
   
    % ------------------------------------------------------------------- %
    % Call readmatrix:
    
        rg_name = strcat('I', num2str(irow_W_YFP), ':', 'I', num2str(irow_W_YFP + nRe - 1));
        W_YFP(its, iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);
        %
        rg_name = strcat('I', num2str(irow_W_mCh), ':', 'I', num2str(irow_W_mCh + nRe - 1));
        W_mCh(its, iBAC, 1:nRe) =  readmatrix(filename, 'Sheet', sheetname, 'Range', rg_name);
        %
        rg_name = strcat('I', num2str(irow_T_YFP), ':', 'I', num2str(irow_T_YFP + nRe - 1));
        T_YFP(its, iBAC, 1:nRe) =  readmatrix(filename, 'Sheet', sheetname, 'Range', rg_name);
        %
        rg_name = strcat('I', num2str(irow_T_mCh), ':', 'I', num2str(irow_T_mCh + nRe - 1));
        T_mCh(its, iBAC, 1:nRe) =  readmatrix(filename, 'Sheet', sheetname, 'Range', rg_name);

        % Actualise indexes of the next BAC concentration:
        irow_W_YFP = irow_W_YFP + nRe;
        irow_W_mCh = irow_W_mCh + nRe;
        irow_T_YFP = irow_T_YFP + nRe;
        irow_T_mCh = irow_T_mCh + nRe;    
    end
    
    % Actualise indexes of the next sampling time:
    irow_W_YFP = irow_ts_W_YFP(its);
    irow_W_mCh = irow_ts_W_mCh(its);
    irow_T_YFP = irow_ts_M_YFP(its);
    irow_T_mCh = irow_ts_M_mCh(its);
end

% ----------------------------------------------------------------------- %
%% (2) Processing data:

% Remove bad replicate of mCherry:
BadRep_W_mCh = [5 3]; % BAC = 75 and Rep = 3 is always = 0;

W_mCh(1:nts, BadRep_W_mCh(1), BadRep_W_mCh(2)) = NaN;

% Join data for the two labels (YFP and mCherry):
W = zeros(nts, nBAC, 2*nRe);
T = W;

W(1:nts, 1:nBAC, 1:nRe)           = W_YFP;
W(1:nts, 1:nBAC, (nRe + 1):2*nRe) = W_mCh;
T(1:nts, 1:nBAC, 1:nRe)           = T_YFP;
T(1:nts, 1:nBAC, (nRe + 1):2*nRe) = T_mCh;

% ----------------------------------------------------------------------- %
%% (3) Calculate average of the replicates:

ave_W_YFP = mean(W_YFP, 3);
ave_W_mCh = nanmean(W_mCh, 3);                                             % Do not account for bad replicate of WT-mCherry;
ave_T_YFP = mean(T_YFP, 3);
ave_T_mCh = mean(T_mCh, 3);

ave_W = nanmean(W, 3);                                                     % Do not account for bad replicate of WT-mCherry;
ave_T = mean(T, 3);

% ----------------------------------------------------------------------- %
%% (4) Calculate survival fractions after 10 min:

% ----------------------------------------------------------------------- %
% (4.a) Differentiating labels:

% Names of the different cases:
strains  = {'W_YFP', 'W_mCh', 'T_YFP', 'T_mCh'};
 
% Names of the survival fractions:
SF_names = {'SF_W_YFP', 'SF_W_mCh', 'SF_T_YFP', 'SF_T_mCh'};

% Number of strains:
nstr     = size(strains, 2);

% Initialise auxiliary survival fraction:
SF       = zeros(nBAC, nRe);

for istr = 1:nstr
    
    % Obtain cell counts for the case:
    eval(sprintf('counts = %s;', strains{istr}))

    % Cell counts at t=0:
    counts_0 = reshape(counts(1, 1:nBAC, 1:nRe), nBAC, nRe);

    % Cell counts at t=10min:
    counts_f = reshape(counts(3, 1:nBAC, 1:nRe), nBAC, nRe);

    for iBAC = 1:nBAC

        % Obtain counts at the current BAC:
        aux_counts_0 = counts_0(iBAC, 1:nRe);
        aux_counts_f = counts_f(iBAC, 1:nRe);

        % Process counts below LDL to avoid problems with log:
        aux_counts_0(aux_counts_0 <= LDL) = 1;
        aux_counts_f(aux_counts_f <= LDL) = 1;

        % Find bad replicates for W-mCherry (NaN) or zero counts at t=0:
        ind_NaN = sort(unique([find(isnan(aux_counts_0)) find(isnan(aux_counts_f)) find(aux_counts_0 == 1)]));

        % Calculate survival fractions:
        SF(iBAC, 1:nRe)   = aux_counts_f./aux_counts_0;
        SF(iBAC, ind_NaN) = NaN;
        
        % Almacenate value for the strain:
        eval(sprintf('%s = SF;', SF_names{istr}))
    end

end

% Save variables without joining labels:
save('Data_TKC.mat', 'nRe', 'ts', 'BAC', 'LDL', strains{:}, SF_names{:})

% ----------------------------------------------------------------------- %
% (4.b) Without differentiating labels:

% Names of the different cases:
strains  = {'W', 'T'};

% Names of the survival fractions:
SF_names = {'SF_W', 'SF_T'};

% Number of cases:
nstr     = size(strains, 2);

% Number of replicates:
nRe      = 2*nRe;

for istr = 1:nstr
    
    % Obtain cell counts for the case:
    eval(sprintf('counts = %s;', strains{istr}))
    
    % Set same cell counts for all replicates at t=0:
    %counts(1, 1:nBAC, 1:nRe) = repmat(reshape(mean(counts(1, 1:nBAC, 1:nRe), 3), 1, nBAC), 1, 1, nRe);
    
    % Cell counts at t=0:
    counts_0 = reshape(counts(1, 1:nBAC, 1:nRe), nBAC, nRe);

    % Cell counts at t=10min:
    counts_f = reshape(counts(3, 1:nBAC, 1:nRe), nBAC, nRe);

    for iBAC = 1:nBAC

        % Obtain counts at the current BAC:
        aux_counts_0 = counts_0(iBAC, 1:nRe);
        aux_counts_f = counts_f(iBAC, 1:nRe);

        % Process counts below LDL:
        aux_counts_0(aux_counts_0 <= LDL) = 1;
        aux_counts_f(aux_counts_f <= LDL) = 1;

        % Find bad replicates for W-mCherry (NaN) or zero counts at t=0:
        ind_NaN = sort(unique([find(isnan(aux_counts_0)) find(isnan(aux_counts_f)) find(aux_counts_0 == 1)]));

        % Calculate survival fractions:
        SF(iBAC, 1:nRe)   = aux_counts_f./aux_counts_0;
        SF(iBAC, ind_NaN) = NaN;
        
        % Almacenate value for the strain:
        eval(sprintf('%s = SF;', SF_names{istr}))
    end

end

% Save data for the strains without differentiating labels:
save('Data_TKC.mat', strains{:}, SF_names{:}, "-append")