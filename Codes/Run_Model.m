%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_Model: Run model predictions for the competition experiment between: 
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
addpath('Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters:

% Join data for YFP and mCherry labels (='_joint') or not (=''):
opts_join = '_joint';

% Number of cycles (can be higher than data):
nc        = 5;

% Duration of the growth periods (min):
t_g       = 24*60;

% Time step for model simulation (min):
ht        = 1/60;

% BAC concentrations of the killing periods (mug/ml):
BAC       = [0 30 40 50].';

% Carrying capacity (CFUS/mL):
K         = 5e8;

% Dilution factor:
D         = 1/100;

% Total initial inocula:
X_0       = 1e6;

% Initial mixing ratios (X_T0/X_W0):
T0_W0     = [1;1e-2;1e-4];

% Culture volume (mL):
V         = 600e-3;

% Extinction concentration (CFUS/mL):
X_e       = 1/V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%% (1) Initialise variables:

% Initialise names of the strains:
if strcmp(opts_join, '_joint')
    strains      = {'W', 'T'};
    ave_SF_names = {'ave_SF_W', 'ave_SF_T'};
    SF_names     = {'SF_W', 'SF_T'};
    CI_SF_names  = {'CI_SF_W', 'CI_SF_T'};
else
    strains      = {'W_YFP', 'W_mCh', 'T_YFP', 'T_mCh'};
    ave_SF_names = {'ave_SF_W_YFP', 'ave_SF_W_mCh', 'ave_SF_T_YFP', 'ave_SF_T_mCh'};
    SF_names     = {'SF_W_YFP', 'SF_W_mCh', 'SF_T_YFP', 'SF_T_mCh'};
    CI_SF_names  = {'CI_SF_W_YFP', 'CI_SF_W_mCh', 'CI_SF_T_YFP', 'CI_SF_T_mCh'};
end

% Simulation times (min):
tsim = (0:ht:nc*t_g).';

% Indexes of the final times of the growth periods:
tcf  = find(ismember(tsim, t_g*(1:nc)));

% Inocula of W and T:
X_W0 = X_0./(1 + T0_W0);
X_T0 = X_0 - X_W0; 

% Problem sizes:
nt     = numel(tsim);                                                      % Number of simulation times:;
nT0_W0 = numel(X_W0);                                                      % Number of initial mixing ratios;
nBAC   = numel(BAC);                                                       % Number of BAC concentrations;
nstr   = size(strains, 2);                                                 % Number of strains;

% Load experimentally determined microbial traits:
load('ProcessData/Data_Traits.mat', 'mu_W', 'mu_T', 'ave_mu_W', 'ave_mu_T',...
     'CI_mu_W', 'CI_mu_T', ave_SF_names{:}, CI_SF_names{:})
 
% Load replicates of the survival fractions:
load('ProcessData/Data_TKC.mat', SF_names{:})

% Get BAC concentrations under interest and add BAC=0:
for istr = 1:nstr
    eval(sprintf('ave_SF_%s = [1;ave_SF_%s];', strains{istr}, strains{istr}))
    eval(sprintf('CI_SF_%s = [1 1;CI_SF_%s];', strains{istr}, strains{istr}))
end

% ----------------------------------------------------------------------- %
%% (2) Simulate model predictions:

% Auxiliary names for the competition cases:
if strcmp(opts_join, '_joint')
    c_names_W = {'W'};
    c_names_T = {'T'};
else
    c_names_W = {'W_YFP', 'W_mCh'};
    c_names_T = {'T_mCh', 'T_YFP'};
end

% Number of competition cases:
n_comp = size(c_names_W, 2);

% ----------------------------------------------------------------------- %
% (2.1) Obtain cell densities of W and T:

% Loop in the competition cases:
for i_comp = 1:n_comp
    
    % Initialise cell densities at the simulation times:
    eval(sprintf('X_%s = zeros(nt, nBAC, nT0_W0);', c_names_W{i_comp}))
    eval(sprintf('X_%s = zeros(nt, nBAC, nT0_W0);', c_names_T{i_comp}))
    
    % Initialise cell densities before disinfection periods:
    eval(sprintf('X_%s_bef = zeros(nc, nBAC, nT0_W0);', c_names_W{i_comp}))
    eval(sprintf('X_%s_bef = zeros(nc, nBAC, nT0_W0);', c_names_T{i_comp}))
    
    % Initialise saturation times:
    eval(sprintf('tsat_%s_%s = zeros(nc, nBAC, nT0_W0);', c_names_W{i_comp}, c_names_T{i_comp}))
    
    % Loop in the different initial mixing ratios:
    for iT0_W0 = 1:nT0_W0
    
        % Obtain initial condition for the current case:
        x_0 = [X_W0(iT0_W0);X_T0(iT0_W0)];

        % Loop in the different BAC concentrations:
        for iBAC = 1:nBAC

            % Obtain survival fractions at the current BAC concentration:
            eval(sprintf('aux_SF_W = ave_SF_%s(iBAC);', c_names_W{i_comp}))
            eval(sprintf('aux_SF_T = ave_SF_%s(iBAC);', c_names_T{i_comp}))
            
            % Call the simulation function:
            input.nc   = nc;
            input.t_g  = t_g;
            input.D    = D;
            input.K    = K;
            input.X_e  = X_e;
            input.mu_W = ave_mu_W;
            input.mu_T = ave_mu_T;
            input.SF_W = aux_SF_W;
            input.SF_T = aux_SF_T;
            input.tsim = tsim;
            input.x_0  = x_0;
            
            output = Sim_Model(input);

            % Save simulated variables:
            eval(sprintf('X_%s(1:nt, iBAC, iT0_W0)       = output.X_W;', c_names_W{i_comp}))
            eval(sprintf('X_%s(1:nt, iBAC, iT0_W0)       = output.X_T;', c_names_T{i_comp}))
            eval(sprintf('X_%s_bef(1:nc, iBAC, iT0_W0)   = output.X_W_bef;', c_names_W{i_comp}))
            eval(sprintf('X_%s_bef(1:nc, iBAC, iT0_W0)   = output.X_T_bef;', c_names_T{i_comp}))
            eval(sprintf('tsat_%s_%s(1:nc, iBAC, iT0_W0) = output.tsat;', c_names_W{i_comp}, c_names_T{i_comp}))
        end
    end
end

% ----------------------------------------------------------------------- %
% (2.2) Obtain extinction cycles:

% Fitness cost of T:
FC = 1 - ave_mu_T/ave_mu_W;
    
% Loop in the competition cases:
for i_comp = 1:n_comp
    
    % Initialise extinction cycles to a high value representing survival
    % for cases where the extinction cycle is not defined:
    eval(sprintf('ce_%s = zeros(nBAC, nT0_W0);', c_names_W{i_comp}))
    eval(sprintf('ce_%s = zeros(nBAC, nT0_W0);', c_names_T{i_comp}))
    
    % Obtain saturation times:
    eval(sprintf('tsat = tsat_%s_%s;', c_names_W{i_comp}, c_names_T{i_comp}))

    % Obtain equilibrium times for the strains at the different BAC conc:
    eval(sprintf('teq_W = - log(D*ave_SF_%s)/ave_mu_W;', c_names_W{i_comp}))
    eval(sprintf('teq_T = - log(D*ave_SF_%s)/ave_mu_T;', c_names_T{i_comp}))

    % Calculate extinction cycles for each BAC concentration:
    for iBAC = 1:nBAC
        
        % Survival advantage of T:
        eval(sprintf('SA = 1 - log(D*ave_SF_%s(iBAC))/log(D*ave_SF_%s(iBAC));', c_names_T{i_comp}, c_names_W{i_comp}))
        
        for iT0_W0 = 1:nT0_W0

            % Obtain saturation time at the first cycle:
            tsat1 = tsat(1, iBAC, iT0_W0);
        
            if ave_mu_W*(tsat1 - teq_W(iBAC)) < log(X_e/X_W0(iT0_W0))      % Extinction in comp. after c=1;
                eval(sprintf('ce_%s(iBAC, iT0_W0) = 1;', c_names_W{i_comp}))
            elseif teq_W(iBAC) > t_g                                       % Extinction of W in isolation;
                if teq_T(iBAC) > t_g || ave_mu_T*(tsat1 - teq_T(iBAC)) < log(X_e/X_T0(iT0_W0)) % Total extinction;
                    eval(sprintf('ce_%s(iBAC, iT0_W0) = ceil((tsat1 - log(X_e/X_W0(iT0_W0))/ave_mu_W - t_g)/(teq_W(iBAC) - t_g));', c_names_W{i_comp}))
                end
            elseif SA - FC > 0
                eval(sprintf('ce_%s(iBAC, iT0_W0) = ceil((log(X_e/X_W0(iT0_W0))/ave_mu_W - log(K/X_T0(iT0_W0))/ave_mu_T + teq_T(iBAC))/(teq_T(iBAC) - teq_W(iBAC)));', c_names_W{i_comp}))
            end
            
            if ave_mu_T*(tsat1 - teq_T(iBAC)) < log(X_e/X_T0(iT0_W0)) % Extinction in comp. after c=1;
                eval(sprintf('ce_%s(iBAC, iT0_W0) = 1;', c_names_T{i_comp}))
            elseif teq_T(iBAC) > t_g % Extinction of W in isolation:
                if teq_W(iBAC) > t_g || ave_mu_W*(tsat1 - teq_W(iBAC)) < log(X_e/X_W0(iT0_W0)) % Total extinction;
                    eval(sprintf('ce_%s(iBAC, iT0_W0) = ceil((tsat1 - log(X_e/X_T0(iT0_W0))/ave_mu_T - t_g)/(teq_T(iBAC) - t_g));', c_names_T{i_comp}))
                end
            elseif SA - FC < 0
                eval(sprintf('ce_%s(iBAC, iT0_W0) = ceil((log(X_e/X_T0(iT0_W0))/ave_mu_T - log(K/X_W0(iT0_W0))/ave_mu_W + teq_W(iBAC))/(teq_W(iBAC) - teq_T(iBAC)));', c_names_T{i_comp}))
            end
        end
    end
end 

% ----------------------------------------------------------------------- %
%% (3) Plot results:

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

% Define colors for different BAC concentrations:
cols     = [183 145 47;
            39 183 222;
            237 106 90;
            129 23 27;
            51 115 87;
            76 57 87;
            227 158 84;
            72 166 167;
            41 115 178]/256;
        
% Define colors for selection and extinction (SE) planes:
col_ext   = [101,115,126]/256;
col_ext_W = col_T;
col_ext_T = col_W;
col_SR_W  = col_W;
col_SR_T  = col_T;
 
% ----------------------------------------------------------------------- %
%% (3.1) Plot traits and confidence intervals:

% Values of x-axis for position of W and T labels:
x  = [1;2];
xL = x(1)-0.5;
xU = x(2)+0.5;

% ----------------------------------------------------------------------- %
% (3.1.a) Plot growth rates:

% Limits of y-axis:
yL = 0;
yU = 1.6;

figure
hold on
set(gcf, 'Color', 'w')
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 25)

% CI for growth rate of W:
fill([x(1)-0.3 x(1)+0.3 x(1)+0.3 x(1)-0.3], 60*[CI_mu_W(1) CI_mu_W(1) CI_mu_W(2) CI_mu_W(2)], col_W, 'EdgeColor', 'k',... 
     'FaceAlpha', 0.6, 'LineWidth', 2) 
 
xaux = x(1)-0.3:0.01:x(1)+0.3;
plot(xaux, 60*ave_mu_W*ones(size(xaux)), 'LineWidth', 2, 'Color', col_W) 

% Plot dots with replicates:
rng(10)
xaux = sqrt(0.005)*randn(size(mu_W)) + x(1);
scatter(xaux, 60*mu_W, 100, col_W, 'filled',  'Marker', 'o', 'MarkerEdgeColor', 'k','LineWidth', 1.5)

% CI for growth rate of T:
fill([x(2)-0.3 x(2)+0.3 x(2)+0.3 x(2)-0.3], 60*[CI_mu_T(1) CI_mu_T(1) CI_mu_T(2) CI_mu_T(2)], col_T, 'EdgeColor', 'k',... 
     'FaceAlpha', 0.6, 'LineWidth', 2) 
xaux = x(2)-0.3:0.01:x(2)+0.3;
plot(xaux, 60*ave_mu_T*ones(size(xaux)), 'LineWidth', 2, 'Color', col_T)  
xaux = sqrt(0.005)*randn(size(mu_T)) + x(2);
scatter(xaux, 60*mu_T, 100, col_T, 'filled',  'Marker', 'o', 'MarkerEdgeColor', 'k','LineWidth', 1.5)

% Plot settings:
xlabel('Strains', 'Interpreter', 'Latex', 'FontSize', 25)
ylabel('Growth rate (1/h)', 'Interpreter', 'Latex', 'FontSize', 25)

xlim([xL xU])
ylim([yL yU])

xticks([x(1) x(2)])
yticks([0 0.4 0.8 1.2 1.6])

xticklabels({'W', 'T'})

ytickformat('%.1f')

ax = gca;

ax.XColor     = 'k';
ax.YColor     = 'k';
ax.LineWidth  = 2;
ax.TickLength = [0.02 0.02]; 

% ----------------------------------------------------------------------- %
% (3.1.b) Plot survival fractions:

% Loop in the competition cases:
for i_comp = 1:n_comp
    
    % Names of the competition cases:
    if n_comp > 1
        fig_name = strcat('W-',c_names_W{i_comp}(3:end),'\&','T-', c_names_T{i_comp}(3:end));
    else
        fig_name = 'W \& T';
    end

    figure
    hold on
    set(gcf, 'Color', 'w')
    set(gca, 'TickLabelInterpreter', 'Latex', 'YScale', 'Log', 'FontSize', 25)
    sgtitle(fig_name, 'Interpreter', 'Latex', 'FontSize', 25)
    for iBAC = 1:(nBAC - 1)
        
        % Obtain survival fractions and CIs:
        eval(sprintf('aux_ave_SF_W = ave_SF_%s(iBAC + 1);', c_names_W{i_comp}))
        eval(sprintf('aux_ave_SF_T = ave_SF_%s(iBAC + 1);', c_names_T{i_comp}))
        eval(sprintf('aux_SF_W     = SF_%s(iBAC, :);', c_names_W{i_comp}))
        eval(sprintf('aux_SF_T     = SF_%s(iBAC, :);', c_names_T{i_comp}))
        eval(sprintf('aux_CI_SF_W  = CI_SF_%s(iBAC + 1, :);', c_names_W{i_comp}))
        eval(sprintf('aux_CI_SF_T  = CI_SF_%s(iBAC + 1, :);', c_names_T{i_comp}))

        % Survival fraction of W:
        fill([x(1)-0.3 x(1)+0.3 x(1)+0.3 x(1)-0.3], [aux_CI_SF_W(1) aux_CI_SF_W(1) aux_CI_SF_W(2) aux_CI_SF_W(2)],...
             cols(iBAC + 1,:), 'EdgeColor', 'k', 'FaceAlpha', 0.6, 'LineWidth', 2) 
        
         % Plot lines with average survival fractions:
        xaux = x(1)-0.3:0.01:x(1)+0.3;
        plot(xaux, aux_ave_SF_W*ones(size(xaux)), 'LineWidth', 2, 'Color',  cols(iBAC + 1,:))

        rng(20)
        xaux = sqrt(0.005)*randn(size(aux_SF_W)) + x(1);
        scatter(xaux, aux_SF_W, 100, cols(iBAC + 1,:), 'filled',  'Marker', 'o', 'MarkerEdgeColor', 'k','LineWidth', 1.5)

        %  Survival fraction of T:
        fill([x(2)-0.3 x(2)+0.3 x(2)+0.3 x(2)-0.3], [aux_CI_SF_T(1) aux_CI_SF_T(1) aux_CI_SF_T(2) aux_CI_SF_T(2)],...
             cols(iBAC + 1,:), 'EdgeColor', 'k', 'FaceAlpha', 0.6, 'LineWidth', 2) 
        xaux = x(2)-0.3:0.01:x(2)+0.3;
        plot(xaux, aux_ave_SF_T*ones(size(xaux)), 'LineWidth', 2, 'Color', cols(iBAC + 1,:)) 
        xaux = sqrt(0.005)*randn(size(aux_SF_T)) + x(2);
        scatter(xaux, aux_SF_T, 100, cols(iBAC + 1,:), 'filled',  'Marker', 'o', 'MarkerEdgeColor', 'k','LineWidth', 1.5)

    end

    % Plot settings:
    xlim([xL xU])
    ylim([5e-10 1])
    
    xticks([x(1) x(2)])
    yticks([1e-9 1e-6 1e-3 1])
    
    xlabel('Strains', 'Interpreter', 'Latex', 'FontSize', 25)
    ylabel('Survival fraction (-)', 'Interpreter', 'Latex', 'FontSize', 25)
    
    if n_comp > 1
    xticklabels({strcat(c_names_W{i_comp}(1),'-',c_names_W{i_comp}(3:end)), strcat(c_names_T{i_comp}(1),'-',c_names_T{i_comp}(3:end))})    
    else
        xticklabels({c_names_W{i_comp}, c_names_T{i_comp}})
    end
    ytickformat('%.1f')

    ax = gca;

    ax.XColor     = 'k';
    ax.YColor     = 'k';
    ax.LineWidth  = 2;
    ax.TickLength = [0.02 0.02];
end

%% (3.b) Model predictions vs flow cytometry data:

% Load flow cytometry (FC) data for the different cases:
load('ProcessData/Data_CompExp_FC_Clustered', 'W_YFP_1_0', 'W_YFP_1_2', 'W_YFP_1_4',...
     'W_mCh_1_0', 'W_mCh_1_2', 'W_mCh_1_4');

% Data sizes:
nc_data    = size(W_mCh_1_0.YFP_freq, 1);
nRe        = size(W_mCh_1_0.YFP_freq, 3);

% Process data for plot:
if strcmp(opts_join, '_joint') > 0
    Data_W = zeros(nc_data, nBAC, nT0_W0, 2*nRe);
    Data_T = zeros(nc_data, nBAC, nT0_W0, 2*nRe);
    
    Data_W(1:nc_data, 1:nBAC, 1, 1:nRe) = W_YFP_1_0.YFP_freq;
    Data_W(1:nc_data, 1:nBAC, 2, 1:nRe) = W_YFP_1_2.YFP_freq;
    Data_W(1:nc_data, 1:nBAC, 3, 1:nRe) = W_YFP_1_4.YFP_freq;
    Data_W(1:nc_data, 1:nBAC, 1, (nRe+1):2*nRe) = W_mCh_1_0.mCh_freq;
    Data_W(1:nc_data, 1:nBAC, 2, (nRe+1):2*nRe) = W_mCh_1_2.mCh_freq;
    Data_W(1:nc_data, 1:nBAC, 3, (nRe+1):2*nRe) = W_mCh_1_4.mCh_freq;
    
    Data_T(1:nc_data, 1:nBAC, 1, 1:nRe) = W_YFP_1_0.mCh_freq;
    Data_T(1:nc_data, 1:nBAC, 2, 1:nRe) = W_YFP_1_2.mCh_freq;
    Data_T(1:nc_data, 1:nBAC, 3, 1:nRe) = W_YFP_1_4.mCh_freq;
    Data_T(1:nc_data, 1:nBAC, 1, (nRe+1):2*nRe) = W_mCh_1_0.YFP_freq;
    Data_T(1:nc_data, 1:nBAC, 2, (nRe+1):2*nRe) = W_mCh_1_2.YFP_freq;
    Data_T(1:nc_data, 1:nBAC, 3, (nRe+1):2*nRe) = W_mCh_1_4.YFP_freq;
else
    Data_W_YFP = zeros(nc_data, nBAC, nT0_W0, nRe);
    Data_W_mCh = zeros(nc_data, nBAC, nT0_W0, nRe);
    Data_T_YFP = zeros(nc_data, nBAC, nT0_W0, nRe);
    Data_T_mCh = zeros(nc_data, nBAC, nT0_W0, nRe);
    
    Data_W_YFP(1:nc_data, 1:nBAC, 1, 1:nRe) = W_YFP_1_0.YFP_freq;
    Data_W_YFP(1:nc_data, 1:nBAC, 2, 1:nRe) = W_YFP_1_2.YFP_freq;
    Data_W_YFP(1:nc_data, 1:nBAC, 3, 1:nRe) = W_YFP_1_4.YFP_freq;
    Data_W_mCh(1:nc_data, 1:nBAC, 1, 1:nRe) = W_mCh_1_0.mCh_freq;
    Data_W_mCh(1:nc_data, 1:nBAC, 2, 1:nRe) = W_mCh_1_2.mCh_freq;
    Data_W_mCh(1:nc_data, 1:nBAC, 3, 1:nRe) = W_mCh_1_4.mCh_freq;
    
    Data_T_mCh(1:nc_data, 1:nBAC, 1, 1:nRe) = W_YFP_1_0.mCh_freq;
    Data_T_mCh(1:nc_data, 1:nBAC, 2, 1:nRe) = W_YFP_1_2.mCh_freq;
    Data_T_mCh(1:nc_data, 1:nBAC, 3, 1:nRe) = W_YFP_1_4.mCh_freq;
    Data_T_YFP(1:nc_data, 1:nBAC, 1, 1:nRe) = W_mCh_1_0.YFP_freq;
    Data_T_YFP(1:nc_data, 1:nBAC, 2, 1:nRe) = W_mCh_1_2.YFP_freq;
    Data_T_YFP(1:nc_data, 1:nBAC, 3, 1:nRe) = W_mCh_1_4.YFP_freq;
end

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

% Plot results:
for i_comp = 1:n_comp
    
    % Obtain cell counts predicted by model:
    eval(sprintf('X_plot_W = X_%s_bef(1:nc, 1:nBAC, 1:nT0_W0);', c_names_W{i_comp}))
    eval(sprintf('X_plot_T = X_%s_bef(1:nc, 1:nBAC, 1:nT0_W0);', c_names_T{i_comp}))
    
    % Obtain data:
    eval(sprintf('FC_Data_W = %s;', sprintf('Data_%s', c_names_W{i_comp})))
    eval(sprintf('FC_Data_T = %s;', sprintf('Data_%s', c_names_T{i_comp})))
 
    % Call plot function:
    if n_comp > 1
        fig_name = strcat('W-',c_names_W{i_comp}(3:end),'\&','T-', c_names_T{i_comp}(3:end));
    else
        fig_name = 'W \& T';
    end
    input.fig_name  = fig_name;
    input.Model_W   = X_plot_W;
    input.Model_T   = X_plot_T;
    input.FC_Data_W = FC_Data_W;
    input.FC_Data_T = FC_Data_T;
    
    % Call to function plotting the results:
    Plot_Model_vs_Data(input)
end

%% (3.c) Model predictions of the subpopulation dynamics:

plot_count = {[1 4 7 10], [2 5 8 11], [3 6 9 12]};

% Plot results:
for i_comp = 1:n_comp
    
    % Names of the competition cases:
    if n_comp > 1
        fig_name = strcat('W-',c_names_W{i_comp}(3:end),'\&','T-', c_names_T{i_comp}(3:end));
    else
        fig_name = 'W \& T';
    end

    % Obtain cell counts predicted by model:
    eval(sprintf('X_plot_W = X_%s(1:nt, 1:nBAC, 1:nT0_W0);', c_names_W{i_comp}))
    eval(sprintf('X_plot_T = X_%s(1:nt, 1:nBAC, 1:nT0_W0);', c_names_T{i_comp}))
    
    fig = figure;
    set(gcf, 'Color', 'w')
    sgtitle(fig_name, 'Interpreter', 'Latex', 'FontSize', 15)

    for iT0_W0 = 1:nT0_W0
        for iBAC = 1:nBAC
            subplot(nBAC, nT0_W0, plot_count{iT0_W0}(iBAC))
            hold on
            set(gca, 'Yscale', 'Log', 'TickLabelInterpreter', 'Latex', 'FontSize', 17)

            plot(tsim, X_plot_W(1:nt, iBAC, iT0_W0), '-', 'Color', col_W, 'LineWidth', 2.5)
            plot(tsim, X_plot_T(1:nt, iBAC, iT0_W0), '-', 'Color', col_T, 'LineWidth', 2.5)

            % ----------------------------------------- %
            % Settings for the subplot:
            if iT0_W0 == 1
                title(sprintf('BAC = %u $(1:10^0)$', BAC(iBAC)), 'Interpreter', 'Latex', 'FontSize', 15)
            elseif iT0_W0 == 2
                title(sprintf('BAC = %u $(1:10^2)$', BAC(iBAC)), 'Interpreter', 'Latex', 'FontSize', 15)
            else
                title(sprintf('BAC = %u $(1:10^4)$', BAC(iBAC)), 'Interpreter', 'Latex', 'FontSize', 15)
            end
            xticks(tsim(tcf))
            xticklabels(1:nc)
            xlim([0 tsim(end)])
            ylim([X_e K])

        end
    end
    
    % Legend:
    if n_comp > 1
        lgd = {strcat('W-',c_names_W{i_comp}(3:end)),strcat('T-', c_names_T{i_comp}(3:end))};
    else
        lgd = {'W', 'T'};
    end
    legend(lgd, 'Location', 'Best', 'Interpreter', 'Latex', 'FontSize', 17)
    
    % Only one xy axes for the figure:
    han = axes(fig,'visible','off'); 

    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';

    Xlim = xlim;           
    Xlb  = mean(Xlim); 
    xlabel(han,{'Cycles'}, 'Interpreter','Latex','FontSize', 17, 'Position',[Xlb -0.07], 'HorizontalAlignment', 'center')
    ylabel(han,{'Cell density (CFU/mL)',''}, 'Interpreter','Latex','FontSize', 17);

end


%% (3.d) Model predictions for the Selection and Extinction (SE) planes:

% Call to the function plotting SE planes:
input.X_0       = X_0;
input.col_ext_W = col_ext_W;
input.col_ext_T = col_ext_T;
input.col_SR_W  = col_SR_W;
input.col_SR_T  = col_SR_T;
input.col_ext   = col_ext;
input.cols      = cols;
input.CI_mu_W   = CI_mu_W;
input.CI_mu_T   = CI_mu_T;

for i_comp = 1:n_comp
    
    % Define inputs for the case:
    eval(sprintf('input.SF_W = ave_SF_%s;', c_names_W{i_comp}))
    eval(sprintf('input.SF_T = ave_SF_%s;', c_names_T{i_comp}))
    eval(sprintf('input.CI_SF_W = CI_SF_%s;', c_names_W{i_comp}))
    eval(sprintf('input.CI_SF_T = CI_SF_%s;', c_names_T{i_comp}))
    
    % Loop in the different initial mixing ratios:
    for iT0_W0 = 1:nT0_W0
        input.T0_W0 = T0_W0(iT0_W0);
        Plot_SE_Plane(input)
        % Names of the competition cases:
        if n_comp > 1
            fig_name = sprintf('Case: %s, Initial mixing ratio = %.4f', strcat('W-',c_names_W{i_comp}(3:end),'\&','T-', c_names_T{i_comp}(3:end)), T0_W0(iT0_W0));
        else
            fig_name = sprintf('Case: %s, Initial mixing ratio = %.4f',  'W \& T', T0_W0(iT0_W0));
        end
        sgtitle(fig_name, 'FontSize', 17, 'Interpreter', 'Latex')
    end
end

%% (3.e) Extinction cycles

% Values of the x-axis:
tplot     = [1 2 3].';

% Value to represent survival of the strains:
surv_line = 0;
for i_comp = 1:n_comp
    eval(sprintf('surv_line = max(max(surv_line, [reshape(ce_%s, [], 1);reshape(ce_%s, [], 1)]));', c_names_W{i_comp}, c_names_T{i_comp}))
end
sur_line  = surv_line + 10;

% Values of the x-axis to plot survival line:
tsurv     = (0:0.1:tplot(end)).';

for i_comp = 1:n_comp
    
    % Obtain extinction cycle of W:
    eval(sprintf('aux_ce_W = ce_%s;', c_names_W{i_comp}));
    
    % Obtain extinction cycle of T:
    eval(sprintf('aux_ce_T = ce_%s;', c_names_T{i_comp}));
       
    figure
    set(gcf, 'Color', 'w', 'Position', [739.0000   96.3333  478.0000  442.6667])
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 25)
    if n_comp > 1
        title(sprintf('Extinction cycle of: %s', strcat(c_names_W{i_comp}(1),'-',c_names_W{i_comp}(3:end))), 'Interpreter', 'Latex', 'FontSize', 25)
    else
        title(sprintf('Extinction cycle of: %s', c_names_W{i_comp}), 'Interpreter', 'Latex', 'FontSize', 25)
    end
    plot(tsurv, surv_line*ones(size(tsurv)), 'k--', 'LineWidth', 2)
    text(tplot(1), surv_line - 3, 'Survival', 'Interpreter', 'Latex', 'FontSize', 20) 
    for iBAC = 1:nBAC
        ind_surv = find(aux_ce_W(iBAC, 1:nT0_W0)==0);
        aux_ce_W(iBAC, ind_surv) = surv_line + iBAC; %#ok<SAGROW>
        plot(tplot, flip(aux_ce_W(iBAC, 1:nT0_W0)), '-o', 'LineWidth', 2, 'Color', cols(iBAC, :),...
             'MarkerFaceColor', cols(iBAC, :), 'MarkerSize', 12)
    end
    xticks(tplot)
    yticks(0:10:(sur_line - 10))
    xticklabels({'$1:10^4$', '$1:10^2$', '$1:1$'})
    xlim([tplot(1)-0.2 tplot(end)+0.2])
    ylim([0 surv_line+10])
    ax = gca;
    ax.XColor    = 'k';
    ax.YColor    = 'k';
    ax.LineWidth = 2;
    ax.TickLength = [0.02 0.02];

    
    figure
    set(gcf, 'Color', 'w', 'Position', [739.0000   96.3333  478.0000  442.6667])
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 25)
    if n_comp > 1
        title(sprintf('Extinction cycle of: %s', strcat(c_names_T{i_comp}(1),'-',c_names_T{i_comp}(3:end))), 'Interpreter', 'Latex', 'FontSize', 25)
    else
        title(sprintf('Extinction cycle of: %s', c_names_T{i_comp}), 'Interpreter', 'Latex', 'FontSize', 25)
    end
    plot(tsurv, surv_line*ones(size(tsurv)), 'k--', 'LineWidth', 2)
    text(tplot(1), surv_line - 3, 'Survival', 'Interpreter', 'Latex', 'FontSize', 20) 
    for iBAC = 1:nBAC
        ind_surv = find(aux_ce_T(iBAC, 1:nT0_W0)==0);
        aux_ce_T(iBAC, ind_surv) = surv_line + iBAC;
        plot(tplot, flip(aux_ce_T(iBAC, 1:nT0_W0)), '-^', 'LineWidth', 2, 'Color', cols(iBAC, :),...
             'MarkerFaceColor', cols(iBAC, :), 'MarkerSize', 12)       
    end
    xticks(tplot)
    yticks(0:10:(sur_line - 10))
    xticklabels({'$1:10^4$', '$1:10^2$', '$1:1$'})
    xlim([tplot(1)-0.2 tplot(end)+0.2])
    ylim([0 surv_line+10])
    ax = gca;
    ax.XColor    = 'k';
    ax.YColor    = 'k';
    ax.LineWidth = 2;
    ax.TickLength = [0.02 0.02];
end

rmpath('Functions')

