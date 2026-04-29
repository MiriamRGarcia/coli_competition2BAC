%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering_GMM: Function implementing clustering (based on Gaussian
% Mixture Model, GMM) on the flow cytometry (FC) data of the competition
% experiment between: 
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
function output = ClusteringGMM(input) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% well_name = Name of the well to perform clustering (e.g., 'A1');
% GMM       = gmdistribution containing the trained Gaussian Mixture Model;
% P_thold   = Interval of probabilities to be member of two clusters;
% Data      = Data of the well to perform clustering;
% cycles    = Array with the cycles to perform clustering;
% plot_res  = Plot clustering results (=1) or not (=0);
% cols      = Colours to plot clusters;
%
% OUTPUT:
% nClus_YFP = Number of events assigned to YFP cluster;
% nClus_mCh = Number of events assigned to mCherry cluster;
% nBack     = Number of events assigned to background noise;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% (1) Initialise variables:

% Call the clustering function:
well_name = input.well;
GMM       = input.GMM;
P_thold   = input.P_thold;
Data      = input.Data;
cycles    = input.cycles;
plot_res  = input.plot_res;
cols      = input.cols;

% Number of cycles:
nc        = numel(cycles);

% Initialise relevant clustering variables for each day:
idClus    = cell(nc, 1);      % Cluster identifiers;
idBoth12  = idClus;           % Identifiers of points in Clusters 1 and 2;
idBoth23  = idClus;           % Identifiers of points in Clusters 2 and 3;
idBoth13  = idClus;           % Identifiers of points in Clusters 1 and 3;
nData     = idClus;           % Number of data points for each cycle;

% ----------------------------------------------------------------------- %
% (2) Perform clustering using the GMM:

if plot_res > 0
    figure
    set(gcf, 'Color', 'w', 'Position', [110 94 1686 833])
    sgtitle(sprintf('Well %s: ', well_name), 'Interpreter', 'Latex', 'FontSize', 15)
end

for ic = 1:nc
    
    % Obtain current day:
    i2c = cycles(ic);
    
    % Obtain data in the log scale:
    Z = log10(1 + [Data.FL1_H{i2c} Data.FL5_H{i2c}]);
    
    % Size of the dataset:
    nData_aux = size(Z, 1);

    % Posterior probabilities of component-member membership:
    P_aux = posterior(GMM, Z);

    % Clasification of the points in clusters according to the mixture model:
    idClus_aux = cluster(GMM, Z);

    % Indexes of the points belonging to more than one cluster:
    idBoth12_aux = find(P_aux(:,1) >= P_thold(1) & P_aux(:,1) <= P_thold(2)); 
    idBoth23_aux = find(P_aux(:,2) >= P_thold(1) & P_aux(:,2) <= P_thold(2)); 
    idBoth13_aux = find(P_aux(:,3) >= P_thold(1) & P_aux(:,3) <= P_thold(2));

    if plot_res > 0
        
        % Plot clustering in the current day:
        subplot(1, nc, ic)
        hold on
        set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
        scatter(Z(idClus_aux == 1, 1), Z(idClus_aux == 1, 2), 20, cols(1, :), 'Marker', 'o')
        scatter(Z(idClus_aux == 2, 1), Z(idClus_aux == 2, 2), 20, cols(2, :), 'Marker', 'o')
        scatter(Z(idClus_aux == 3, 1), Z(idClus_aux == 3, 2), 20, cols(3, :), 'Marker', 'x')
 
        % Points belonging to cluster boundaries:
        plot(Z(idBoth12_aux,1), Z(idBoth12_aux,2), 'ko', 'MarkerSize', 5)
        plot(Z(idBoth23_aux,1), Z(idBoth23_aux,2), 'ko', 'MarkerSize', 5)
        plot(Z(idBoth13_aux,1), Z(idBoth13_aux,2), 'ko', 'MarkerSize', 5)

        xlabel('FITC-H', 'Interpreter', 'Latex', 'FontSize', 13)
        ylabel('PE-H', 'Interpreter', 'Latex', 'FontSize', 13)
        xlim([2 4.5])
        ylim([2 4.5])
        xticks([1.5 2 2.5 3 3.5 4 4.5])
        legend({'mCherry', 'YFP', 'Background', 'Boundary'}, 'Location', 'NorthWest', 'Interpreter', 'Latex', 'FontSize', 15, 'EdgeColor', 'None')
        hold off
    end
    
    % Save relevant clustering variables:
    idClus{ic}   = idClus_aux;
    idBoth12{ic} = idBoth12_aux;
    idBoth23{ic} = idBoth23_aux;
    idBoth13{ic} = idBoth13_aux;
    nData{ic}    = nData_aux;
end

% ----------------------------------------------------------------------- %
% (3) Obtain the composition of the population:

if plot_res > 0
    figure
    set(gcf, 'Color', 'w', 'Position', [110 94 1686 833])
    sgtitle(sprintf('Well %s: ', well_name), 'Interpreter', 'Latex', 'FontSize', 15)
end

% Number of elements of the clusters for the different cycles:
nClus_YFP = zeros(nc, 1);
nClus_mCh = nClus_YFP;
nBack     = nClus_YFP;

for ic = 1:nc

    % Obtain the points in the YFP/mCherry/Background clusters:
    idClus_mCh  = idClus{ic} == 1;
    idClus_YFP  = idClus{ic} == 3;
    idClus_Back = idClus{ic} == 2;

    % Remove points belonging to more than one cluster:
    ii_Both             = sort(unique([idBoth12{ic};idBoth23{ic};idBoth13{ic}]));
    idClus_mCh(ii_Both)  = 0;
    idClus_YFP(ii_Both)  = 0;
    idClus_Back(ii_Both) = 0;

    % Size of the clusters:
    nClus_YFP(ic) = numel(find(idClus_YFP));
    nClus_mCh(ic) = numel(find(idClus_mCh)); 
    nBack(ic)     = numel(find(idClus_Back));

    % Calculate the percentage of each subpopulation:
    FracClus_YFP = nClus_YFP(ic)/(nClus_mCh(ic) + nClus_YFP(ic));
    FracClus_mCh = nClus_mCh(ic)/(nClus_mCh(ic) + nClus_YFP(ic));
    
    % Obtain data in the current cycle:
    i2c = cycles(ic);
    
    % Obtain data in the log scale:
    Z = log10(1 + [Data.FL1_H{i2c} Data.FL5_H{i2c}]);
    
    if plot_res > 0
        % Plot clustering and fractions:
        subplot(1, nc, ic)
        hold on
        set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 13)
        scatter(Z(idClus_mCh, 1), Z(idClus_mCh, 2), 20, cols(1, :))
        scatter(Z(idClus_YFP, 1), Z(idClus_YFP, 2), 20, cols(3, :))

        xlabel('FITC-H', 'Interpreter', 'Latex', 'FontSize', 13)
        ylabel('PE-H', 'Interpreter', 'Latex', 'FontSize', 13)
        xlim([2 4.5])
        ylim([2 4.5])
        xticks([1.5 2 2.5 3 3.5 4 4.5])
        auxlgd = strcat(sprintf('mCherry $=%.2f', 100*FracClus_mCh), '\%$');
        lgd    = {auxlgd, strcat(sprintf('YFP $=%.2f', 100*FracClus_YFP), '\%$')};
        legend(lgd, 'FontSize', 15, 'EdgeColor', 'None', 'Interpreter','latex')    
        
    end

end

output.nClus_YFP = nClus_YFP;
output.nClus_mCh = nClus_mCh;
output.nBack     = nBack;

if plot_res > 0
    pause
    fprintf('\n >> Press any button to continue.')
    close all
end

end