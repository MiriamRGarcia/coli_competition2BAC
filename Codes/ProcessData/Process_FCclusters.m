%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process_FCclusters: Obtain frequencies of YFP and mCherry after
% clustering the flow cytometry (FC) data of the competition experiment
% between: 
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
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% Minimum number of events per microliter to calculate statistics:
min_nEv   = 100;

% Plot results of the clustering algorithm:
plot_res  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Initialise variables:

% Names of the strains:
strains   = {'W_mCh_1_0', 'W_mCh_1_2', 'W_mCh_1_4',...
             'W_YFP_1_0', 'W_YFP_1_2', 'W_YFP_1_4',...
             'W_mCh', 'T_mCh'};
         
% Problem sizes:
nstr   = size(strains, 2);
nRe    = 3;
nBAC   = 4;
nc     = 5;

% Array with cycles:
cycles = 1:nc;

% Letters of all the wells:
let_names = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};

% Row numbers on the 96-well plate for all the cases:
row       = [1 2 3;7 8 9;1 2 3;4 5 6;10 11 12;4 5 6;7 8 9;10 11 12];

% Identify if the cases are in the top (=0) or bottom (=1) of the plate:
id_row    = [0;0;1;0;0;1;1;1];

% Sample volume recorded with FC each cycle:
rec_mL    = [5;5;5;5;3];

% Minimum number of events for each cycle:
min_nEv   = rec_mL*min_nEv;

% ----------------------------------------------------------------------- %
%% (2) Run clustering for each case:

% Load trained Gaussian Mixture Model (GMM):
load('Data_CompExp_TrainedGMM.mat')

% Colors for clusters:
cols = [229 127 127;
        38 91 143;
        207 179 65]/256;
    
for istr = 1:nstr
    
    % Initialise variables:
    eval(sprintf('%s.YFP_nData = zeros(nc, nBAC, nRe);', strains{istr}))
    eval(sprintf('%s.mCh_nData = zeros(nc, nBAC, nRe);', strains{istr}))
    eval(sprintf('%s.YFP_freq  = zeros(nc, nBAC, nRe);', strains{istr}))
    eval(sprintf('%s.mCh_freq  = zeros(nc, nBAC, nRe);', strains{istr}))
    
    % Detect if the case is in the top of the plate:
    if id_row(istr) < 1 % The case is in the top of the plate;
        
        for iRe = 1:nRe

            for iBAC = 1:nBAC
                
                % Name of the well to perform clustering:
                well_name = strcat(let_names{iBAC}, num2str(row(istr, iRe)));
                
                % Load FC data of the well from .mat file:
                load('Data_CompExp_FC.mat', well_name)
                
                % Create temporal variable with the well data:
                eval(sprintf('Data = %s;', well_name))
                
                % Call the clustering function (classify data into
                % one of the trained GMM clusters):
                input.strain   = strains{istr};
                input.well     = well_name;
                input.GMM      = GMM;
                input.P_thold  = P_thold;
                input.Data     = Data;
                input.cycles   = cycles;
                input.plot_res = plot_res;
                input.cols     = cols;
                
                output         = ClusteringGMM(input);
                
                % Obtain number of events assigned to each cluster:
                nClus_YFP = output.nClus_YFP;
                nClus_mCh = output.nClus_mCh;
                
                % Save variables:       
                nClus = nClus_YFP + nClus_mCh;
                eval(sprintf('%s.YFP_nData(1:nc, %u, %u) = nClus_YFP;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.mCh_nData(1:nc, %u, %u) = nClus_mCh;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.YFP_freq(1:nc, %u, %u) = nClus_YFP./nClus;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.mCh_freq(1:nc, %u, %u) = nClus_mCh./nClus;', strains{istr}, iBAC, iRe))
                
                % Detect cycles with cero events:
                ext_flag = find(nClus - min_nEv < 0);
                eval(sprintf('%s.YFP_freq(ext_flag, %u, %u) = NaN;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.mCh_freq(ext_flag, %u, %u) = NaN;', strains{istr}, iBAC, iRe))             
            end
        end
    else           % The case is in the bottom of the plate;
        for iRe = 1:nRe
            
            for iBAC = 1:nBAC
                
                % Name of the well for clustering:
                well_name = strcat(let_names{nBAC + iBAC}, num2str(row(istr, iRe)));
                
                % Load FC data of the well from .mat file:
                load('Data_CompExp_FC', well_name)
                
                % Create temporal variable with the well data:
                eval(sprintf('Data = %s;', well_name))
                
                % Call the clustering function:
                input.strain   = strains{istr};
                input.well     = well_name;
                input.GMM      = GMM;
                input.Data     = Data;
                input.cycles   = cycles;
                input.plot_res = plot_res;
                input.cols     = cols;
                
                output         = ClusteringGMM(input);
                
                % Obtain number of events assigned to each cluster:
                nClus_YFP = output.nClus_YFP;
                nClus_mCh = output.nClus_mCh;
                
                % Save variables:     
                nClus = nClus_YFP + nClus_mCh;
                eval(sprintf('%s.YFP_nData(1:nc, %u, %u) = nClus_YFP;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.mCh_nData(1:nc, %u, %u) = nClus_mCh;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.YFP_freq(1:nc, %u, %u) = nClus_YFP./nClus;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.mCh_freq(1:nc, %u, %u) = nClus_mCh./nClus;', strains{istr}, iBAC, iRe))
                
                % Detect cycles with cero events:
                ext_flag = find(nClus - min_nEv < 0);
                eval(sprintf('%s.YFP_freq(ext_flag, %u, %u) = NaN;', strains{istr}, iBAC, iRe))
                eval(sprintf('%s.mCh_freq(ext_flag, %u, %u) = NaN;', strains{istr}, iBAC, iRe))        
            end
        end
    end
end


% ----------------------------------------------------------------------- %
%% (3) Calculate average of the replicates:

for istr = 1:nstr
    
    eval(sprintf('aux = %s;', strains{istr}))
    
    % Initialise averages:
    aux_ave_YFP = zeros(nc, nBAC);
    aux_ave_mCh = zeros(nc, nBAC);
            
    for ic = 1:nc
        for iBAC = 1:nBAC 
            ngoodRe = 0;
            for iRe = 1:nRe
                if ~isnan(aux.YFP_freq(ic, iBAC, iRe))
                    aux_ave_YFP(ic, iBAC) = aux_ave_YFP(ic, iBAC) + aux.YFP_freq(ic, iBAC, iRe);
                    ngoodRe = ngoodRe + 1;
                end
            end
            aux_ave_YFP(ic, iBAC) = aux_ave_YFP(ic, iBAC)/ngoodRe;
            %
            ngoodRe = 0;
            for iRe = 1:nRe
                if ~isnan(aux.mCh_freq(ic, iBAC, iRe))
                    aux_ave_mCh(ic, iBAC) = aux_ave_mCh(ic, iBAC) + aux.mCh_freq(ic, iBAC, iRe);
                    ngoodRe = ngoodRe + 1;
                end
            end
            aux_ave_mCh(ic, iBAC) = aux_ave_mCh(ic, iBAC)/ngoodRe;
            
            eval(sprintf('%s.YFP_freq_ave = aux_ave_YFP;', strains{istr}))
            eval(sprintf('%s.mCh_freq_ave = aux_ave_mCh;', strains{istr}))
        end
    end
end
        
% ----------------------------------------------------------------------- %
% Save results:
save('Data_CompExp_FC_Clustered.mat', strains{:})