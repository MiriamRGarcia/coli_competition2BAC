%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read_CompExp: Read data from the Competition Experiment (CompExp) 
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
% Cell count data were collected before each disinfection round.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% Number of cycles for the competition experiment:
nc          = 4;

% Number of replicates of the competition experiment (without joining
% labels):
nRe         = 3;

% Exposure time for the competition experiment (min):
tk          = 10;

% BAC concentrations for the competition experiment (mug/ml):
BAC         = [0 30 40 50].';

% Plot results of processing FC data (=1) or not (=0):
plot_res    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Read cell count data from Excel:

% Name of the Excel file:
filename = '../../Data/CellCounts_CompExp_WvsT_YFP_mCherry_BAC_M9.xlsx';

% Initialise cell counts, where:
% 1st element: Time;
% 2nd element: BAC concentration;
% 3th element: Replicate;
W_mCh_1_0 = cell(nc, 1);
W_mCh_1_2 = W_mCh_1_0;
W_mCh_1_4 = W_mCh_1_0;
W_YFP_1_0 = W_mCh_1_0;
W_YFP_1_2 = W_mCh_1_0;
W_YFP_1_4 = W_mCh_1_0;
W_mCh     = W_mCh_1_0;
T_mCh     = W_mCh_1_0;

% Problem sizes:
nBAC = numel(BAC);

% Auxiliary variables:
W_mCh_1_0_aux = zeros(nBAC, nRe);
W_mCh_1_2_aux = W_mCh_1_0_aux;
W_mCh_1_4_aux = W_mCh_1_0_aux;
W_YFP_1_0_aux = W_mCh_1_0_aux;
W_YFP_1_2_aux = W_mCh_1_0_aux;
W_YFP_1_4_aux = W_mCh_1_0_aux;
W_mCh_aux     = W_mCh_1_0_aux;
T_mCh_aux     = W_mCh_1_0_aux;

ave_W_mCh_1_0 = W_mCh_1_0;
ave_W_mCh_1_2 = W_mCh_1_0;
ave_W_mCh_1_4 = W_mCh_1_0;
ave_W_YFP_1_0 = W_mCh_1_0;
ave_W_YFP_1_2 = W_mCh_1_0;
ave_W_YFP_1_4 = W_mCh_1_0;
ave_W_mCh     = W_mCh_1_0;
ave_T_mCh     = W_mCh_1_0;

% Loop in the number of cycles:
for ic = 1:nc
    sheetname = sprintf('cycle%u', ic);
    
    % Initialise row counters:
    irow_W_mCh_1_0 = 2;
    irow_W_mCh_1_2 = 26;
    irow_W_mCh_1_4 = 50;
    irow_W_YFP_1_0 = 14;
    irow_W_YFP_1_2 = 38;
    irow_W_YFP_1_4 = 62;
    irow_W_mCh     = 74;
    irow_T_mCh     = 86;

    % Loop in the BAC concentrations:
    for iBAC = 1:nBAC

        % Define row range for sampling time itS and BAC iBAC:
        rg_name = strcat('M', num2str(irow_W_mCh_1_0), ':', 'M', num2str(irow_W_mCh_1_0 + nRe - 1));

        % Call readmatrix:
        W_mCh_1_0_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);

        rg_name = strcat('M', num2str(irow_W_mCh_1_2), ':', 'M', num2str(irow_W_mCh_1_2 + nRe - 1));
        W_mCh_1_2_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);

        rg_name = strcat('M', num2str(irow_W_mCh_1_4), ':', 'M', num2str(irow_W_mCh_1_4 + nRe - 1));
        W_mCh_1_4_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);

        rg_name = strcat('M', num2str(irow_W_YFP_1_0), ':', 'M', num2str(irow_W_YFP_1_0 + nRe - 1));
        W_YFP_1_0_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);

        rg_name = strcat('M', num2str(irow_W_YFP_1_2), ':', 'M', num2str(irow_W_YFP_1_2 + nRe - 1));
        W_YFP_1_2_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);            

        rg_name = strcat('M', num2str(irow_W_YFP_1_4), ':', 'M', num2str(irow_W_YFP_1_4 + nRe - 1));
        W_YFP_1_4_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);

        rg_name = strcat('M', num2str(irow_W_mCh), ':', 'M', num2str(irow_W_mCh + nRe - 1));
        W_mCh_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);            

        rg_name = strcat('M', num2str(irow_T_mCh), ':', 'M', num2str(irow_T_mCh + nRe - 1));
        T_mCh_aux(iBAC, 1:nRe) = readmatrix(filename,'Sheet', sheetname, 'Range', rg_name);

        % Actualise indexes of the next BAC concentration:
        irow_W_mCh_1_0 = irow_W_mCh_1_0 + nRe;
        irow_W_mCh_1_2 = irow_W_mCh_1_2 + nRe;
        irow_W_mCh_1_4 = irow_W_mCh_1_4 + nRe;
        irow_W_YFP_1_0 = irow_W_YFP_1_0 + nRe;
        irow_W_YFP_1_2 = irow_W_YFP_1_2 + nRe;
        irow_W_YFP_1_4 = irow_W_YFP_1_4 + nRe;
        irow_W_mCh     = irow_W_mCh + nRe;
        irow_T_mCh     = irow_T_mCh + nRe;

    end
    
    % Almacenate data for the different days:
    W_mCh_1_0{ic} = W_mCh_1_0_aux;
    W_mCh_1_2{ic} = W_mCh_1_2_aux;
    W_mCh_1_4{ic} = W_mCh_1_4_aux;
    W_YFP_1_0{ic} = W_YFP_1_0_aux;
    W_YFP_1_2{ic} = W_YFP_1_2_aux;
    W_YFP_1_4{ic} = W_YFP_1_4_aux;
    W_mCh{ic}     = W_mCh_aux;
    T_mCh{ic}     = T_mCh_aux;
       
    % Calculate average of the replicates:
    ave_W_mCh_1_0{ic} = mean(W_mCh_1_0_aux, 2);
    ave_W_mCh_1_2{ic} = mean(W_mCh_1_2_aux, 2);
    ave_W_mCh_1_4{ic} = mean(W_mCh_1_4_aux, 2);
    ave_W_YFP_1_0{ic} = mean(W_YFP_1_0_aux, 2);
    ave_W_YFP_1_2{ic} = mean(W_YFP_1_2_aux, 2);
    ave_W_YFP_1_4{ic} = mean(W_YFP_1_4_aux, 2);
    ave_W_mCh{ic}     = mean(W_mCh_aux, 2);
    ave_T_mCh{ic}     = mean(T_mCh_aux, 2);
    
        % Restore indexes to the original value:
        irow_W_mCh_1_0 = irow_W_mCh_1_0 + nRe;
        irow_W_mCh_1_2 = irow_W_mCh_1_2 + nRe;
        irow_W_mCh_1_4 = irow_W_mCh_1_4 + nRe;
        irow_W_YFP_1_0 = irow_W_YFP_1_0 + nRe;
        irow_W_YFP_1_2 = irow_W_YFP_1_2 + nRe;
        irow_W_YFP_1_4 = irow_W_YFP_1_4 + nRe;
        irow_W_mCh     = irow_W_mCh + nRe;
        irow_T_mCh     = irow_T_mCh + nRe;
end

% Save cell count data:
save('Data_CompExp_CellCounts.mat', 'BAC', 'nRe', 'nc', 'W_mCh', 'T_mCh', ...
     'W_mCh_1_0', 'W_mCh_1_2', 'W_mCh_1_4', 'W_YFP_1_0', 'W_YFP_1_2', 'W_YFP_1_4', ...
     'ave_W_mCh_1_0', 'ave_W_mCh_1_2', 'ave_W_mCh_1_4', 'ave_W_YFP_1_0', 'ave_W_YFP_1_2', ...
     'ave_W_YFP_1_4', 'ave_W_mCh', 'ave_T_mCh');

% -----------------------------------------------------------------------
%% (2) Read FC data:

clearvars -except nc plot_res

% Initialise information of the wells:
let_well = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
num_well = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
nl       = size(let_well, 2);
nw       = size(num_well, 2);

for il = 1:nl
    for iw = 1:nw
        
        % Name of the well to read:
        well_name = strcat(let_well{il}, num_well{iw});
        
        for ic = 1:(nc + 1)

            % Name of the file containing the data for the well:
            filename = sprintf('../../Data/FlowCytometry_CompExp_WvsT_YFP_mCherry_BAC_M9/Cycle_%u/01-Well-%s.fcs', ic, well_name);
            
            % Process FC data of the current well (remove noise and doublets):
            auxData = Process_FCdata(filename, well_name, ic, plot_res);  %#ok<NASGU>

            % Almacenate relevant variables for the well:
            eval(sprintf('%s.FSC_H{%u} = auxData.FSC_H;', well_name, ic))
            eval(sprintf('%s.FSC_A{%u} = auxData.FSC_A;', well_name, ic))
            eval(sprintf('%s.SSC_H{%u} = auxData.SSC_H;', well_name, ic))
            eval(sprintf('%s.SSC_A{%u} = auxData.SSC_A;', well_name, ic))
            eval(sprintf('%s.FL1_H{%u} = auxData.FL1_H;', well_name, ic))
            eval(sprintf('%s.FL1_A{%u} = auxData.FL1_A;', well_name, ic))
            eval(sprintf('%s.FL5_H{%u} = auxData.FL5_H;', well_name, ic))
            eval(sprintf('%s.FL5_A{%u} = auxData.FL5_A;', well_name, ic))
            close
        end
    end
end

% Save FC data:
clear well_name num_well let_well il iw ic nc nl nw filename auxData plot_res
  
save('Data_CompExp_FC.mat')
