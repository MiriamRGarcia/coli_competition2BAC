%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot_Model_vs_Data: Function plotting the results of model-based simulation
% of the competition experiment, and comparing with flow cytometry data
% Strains:
% Wild-type E. coli YFP (W-YFP) vs Tolerant E. coli mCherry (T-mCh),
% Wild-type E. coli mCherry (W-mCh) vs Tolerant E. coli YFP (T-YFP).
% The competing strains were exposed to sucessive cycles alternating 
% growth periods and disinfection periods with benzalkonium chloride. 
% After each disinfection step, the cultures were diluted into 
% fresh medium to start the next treatment cycle (4 complete cycles).
% Growth periods lasted 24 hours.
% Disinfection periods lasted 10 minutes.
% Flow cytometry (FC) data was collected before each disinfection round.
% Cell count data were collected before and after each disinfection round.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_Model_vs_Data(input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% (1) Initialise variables:

% Obtain input variables:
fig_name   = input.fig_name;
tplot      = input.tplot;
tplot_data = input.tplot_data;
tlabs      = input.tlabs;
BAC        = input.BAC;
Model_W    = input.Model_W;
Model_T    = input.Model_T;
FC_Data_W  = input.FC_Data_W;
FC_Data_T  = input.FC_Data_T;
X_e        = input.X_e;
col_W      = input.col_W;
col_T      = input.col_T;
col_Re_W   = input.col_Re_W;
col_Re_T   = input.col_Re_T;

% Problem sizes:
nc      = size(Model_W, 1);
nc_data = numel(tplot_data);
nBAC    = size(Model_W, 2);
nT0_W0  = size(Model_W, 3);
nRe     = size(FC_Data_W, 4);
          
% Define counters for subplots:
plot_count = {[1 4 7 10], [2 5 8 11], [3 6 9 12]};
          
% ----------------------------------------------------------------------- %
% (2) Plot results:

fig = figure;
set(gcf, 'Color', 'w')
sgtitle(fig_name, 'Interpreter', 'Latex', 'FontSize', 15)

for iT0_W0 = 1:nT0_W0
    for iBAC = 1:nBAC
        
        % Obtain model predictions:
        Model_W_aux = Model_W(1:nc, iBAC, iT0_W0);
        Model_T_aux = Model_T(1:nc, iBAC, iT0_W0);
        
        % Average total cell counts before killing:
        Model_aux   = Model_W_aux + Model_T_aux; 
        
        % Plot results:
        subplot(nBAC, nT0_W0, plot_count{iT0_W0}(iBAC));
        hold on
        set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 17)
        
        % -------------------------------------------- %
        % (2.a) Plot model predictions:

        % Detect total population extinction:
        indplot = find(max(Model_W_aux, Model_T_aux) < X_e);

        Model_W_aux          = 100*Model_W_aux./Model_aux;
        Model_W_aux(indplot) = 0;

        Model_T_aux = 100*Model_T_aux./Model_aux;
        Model_T_aux(indplot) = 0;
        
        plot(tplot, Model_W_aux, '-', 'Color', col_W, 'LineWidth', 2.5)
        plot(tplot, Model_T_aux, '-', 'Color', col_T, 'LineWidth', 2.5)
        
        % --------------------------------------------------------------- %
        % (2.b) Plot FC data:
        
        % --------------------------------------------------------------- %
        % Fill variability range between replicates:
        Data_aux  = reshape(100*FC_Data_W(1:nc_data, iBAC, iT0_W0, 1:nRe), nc_data, nRe);
        
        % Remove wrong replicate:
        if iBAC > 1
            Data_aux(isnan(Data_aux)) = 0;
        end

        lowlim = min(Data_aux, [], 2);
        uplim  = max(Data_aux, [], 2);
        
        taux = tplot_data;

        coord_up  = [taux, uplim];
        coord_low = [taux, lowlim];
    
        coord_combine = [coord_up;flipud(coord_low)];
    
        fill(coord_combine(:,1), coord_combine(:,2), col_W, 'EdgeColor', 'None',...
             'FaceAlpha', 0.2, 'HandleVisibility', 'off')
        

        Data_aux = reshape(100*FC_Data_T(1:nc_data, iBAC, iT0_W0, 1:nRe), nc_data, nRe);
        
        % Remove wrong replicate:
        if iT0_W0 > 1 || iBAC > 1
             Data_aux(isnan(Data_aux)) = 0;
        end
        
        lowlim = min(Data_aux, [], 2);
        uplim  = max(Data_aux, [], 2);
      
        taux = tplot_data;
        
        coord_up  = [taux, uplim];
        coord_low = [taux, lowlim];
    
        coord_combine = [coord_up;flipud(coord_low)];
    
        fill(coord_combine(:,1), coord_combine(:,2), col_T, 'EdgeColor', 'None',...
             'FaceAlpha', 0.2, 'HandleVisibility', 'off')

        % --------------------------------------------------------------- %
        % Plot data replicates:
        
        for iRe = 1:nRe
            scatter(tplot_data, 100*FC_Data_W(1:nc_data, iBAC, iT0_W0, iRe), 50, col_Re_W(iRe,:), 'filled')
            scatter(tplot_data, 100*FC_Data_T(1:nc_data, iBAC, iT0_W0, iRe), 50, col_Re_T(iRe,:), 'filled', 'Marker', '^')
        end
               
        % Plot asteric for extinct replicates:
        for ic = 1:nc_data
            n_zcounts = numel(find(Data_aux(ic, 1:nRe)==0));
            if n_zcounts > 0
                n_ztext = strcat('*', num2str(n_zcounts));
                text(tplot_data(ic) + 0.1, 8, n_ztext, 'FontSize', 13, 'Interpreter', 'Latex')
            end
        end
        
        hold off

        % --------------------------------------------------------------- %
        % (2.c) Settings for the subplot:

        xticks(tplot)
        xticklabels(tlabs)
        yticks([0 25 50 75 100])
        %xlim([0 nc+1])
        ylim([0 105])
        
        ax = gca;
        
        ax.XColor     = 'k';
        ax.YColor     = 'k';
        ax.LineWidth  = 1.2;
        ax.TickLength = [0.02 0.02];
        
        if iT0_W0 == 1
            title(sprintf('BAC = %u $(1:10^0)$', BAC(iBAC)), 'Interpreter', 'Latex', 'FontSize', 15)
        elseif iT0_W0 == 2
            title(sprintf('BAC = %u $(1:10^2)$', BAC(iBAC)), 'Interpreter', 'Latex', 'FontSize', 15)
        else
            title(sprintf('BAC = %u $(1:10^4)$', BAC(iBAC)), 'Interpreter', 'Latex', 'FontSize', 15)
        end
    end
end

% Only one xy axes for the figure:
han = axes(fig,'visible','off'); 

han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

Xlim = xlim;           
Xlb  = mean(Xlim); 
xlabel(han,{'Cycles'}, 'Interpreter','Latex','FontSize', 17, 'Position',[Xlb -0.07], 'HorizontalAlignment', 'center')
ylabel(han,{'Subopulation frequency (\%)',''}, 'Interpreter','Latex','FontSize', 17);


end