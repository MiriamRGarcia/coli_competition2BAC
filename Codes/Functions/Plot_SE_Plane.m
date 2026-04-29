%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot_SE_Plane: Plot selection and extinction (SE) planes of the
% competition experiment between strains: 
% Wild-type E. coli YFP (W-YFP) vs Tolerant E. coli mCherry (T-mCh),
% Wild-type E. coli mCherry (W-mCh) vs Tolerant E. coli YFP (T-YFP).
% Competition cases W-YFP vs T-mCh and W-mCh vs T-YFP treated together.
% The competing strains were exposed to sucessive cycles alternating 
% growth periods with periods of disinfection with benzalkonium chloride. 
% After disinfection, the cultures were diluted into fresh medium to 
% start the next cycle.
% Growth periods lasted 24 hours.
% Disinfection periods lasted 10 minutes.
% Flow cytometry (FC) data was collected before each disinfection round.
% Cell count data were collected before and after each disinfection round.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_SE_Plane(input) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% (1) Initialise variables:

% Setup parameters:
t_g         = input.t_g;
X_0         = input.X_0;
T0_W0       = input.T0_W0;
K           = input.K;
D           = input.D;
X_e         = input.X_e;

% Microbial traits:
mu_W       = input.mu_W;
mu_T       = input.mu_T;
SF_W       = input.SF_W;
SF_T       = input.SF_T;

% Confidence intervals for the microbial traits:
CI_mu_W    = input.CI_mu_W;
CI_mu_T    = input.CI_mu_T;
CI_SF_W    = input.CI_SF_W;
CI_SF_T    = input.CI_SF_T;

% Tranformation of SF to log-scale:
logSF_W    = log(SF_W);
logSF_T    = log(SF_T);
CI_logSF_W = log(CI_SF_W);
CI_logSF_T = log(CI_SF_T);

% Colours:
cols       = input.cols;
col_ext_W  = input.col_ext_W;
col_ext_T  = input.col_ext_T;
col_SR_W   = input.col_SR_W;
col_SR_T   = input.col_SR_T;
col_ext    = input.col_ext;

% Number of different experiments:
nBAC = numel(SF_W);

% BAC concentrations:
BAC  = input.BAC;

% Initial concentration (CFU/mL) of each strain:
X_W0 = X_0/(1 + T0_W0);
X_T0 = X_0 - X_W0; 

% ----------------------------------------------------------------------- %
% (2) First plot (selection plane):
    
% Calculate fitness cost of T relative to W:
FC = 1 - mu_T/mu_W;

% Bounds on fitness cost if provided:
try
    min_FC = 1 - CI_mu_T(2)/CI_mu_W(1);
    max_FC = 1 - CI_mu_T(1)/CI_mu_W(2);
catch
    min_FC = [];
    max_FC = [];
end

figure;
hold on
set(gcf, 'Color', 'w')
set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 25)

% Plot equilibrium line:
x = linspace(0, 1, 100);
plot(x, x, 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
 
% Plot text indicating selection or not:
text(0.05, 0.9, 'T is selected', 'FontSize', 15, 'Interpreter', 'Latex')
text(0.75, 0.1, 'W is selected', 'FontSize', 15, 'Interpreter', 'Latex')
text(0.45, 0.5, 'Equilibrium', 'FontSize', 15, 'Interpreter','latex', 'Rotation', 40)

% Plot the boxes for each concentration:
lgd = cell(nBAC, 1);
for iBAC = 1:nBAC

    % Calculate survival advantage of T over W:
    SA = 1 - (logSF_T(iBAC) + log(D))/(logSF_W(iBAC) + log(D));

    % Bounds on survival advantage if provided:
    try
        min_SA = 1 - (CI_logSF_T(iBAC, 1) + log(D))/(CI_logSF_W(iBAC, 2) + log(D));
        max_SA = 1 - (CI_logSF_T(iBAC, 2) + log(D))/(CI_logSF_W(iBAC, 1) + log(D));
    catch
        min_SA = [];
        max_SA = [];
    end 

    % Plot box with bounds on FC and SA:
    aux_x = [min_FC max_FC max_FC min_FC];
    aux_y = [min_SA min_SA max_SA max_SA];

    if numel(aux_x) == numel(aux_y) && numel(aux_x) == 4 % Plot box;
        fill(aux_x, aux_y, cols(iBAC, :), 'EdgeColor', cols(iBAC, :), 'FaceAlpha', 0.4,...
            'LineWidth', 2, 'HandleVisibility', 'off')

    elseif numel(aux_x) == 4 % Plot line indicating interval for FC;
        aux_x = linspace(min_FC, max_FC, 100);
        plot(aux_x, SA*ones(size(aux_x)), 'LineWidth', 2, 'Color', cols(iBAC, :), 'HandleVisibility', 'off')

    elseif numel(aux_y) == 4 % Plot line indicating interval for SA;
        aux_y = linspace(min_SA, max_SA, 100);
        plot(FC*ones(size(aux_y)), aux_y, 'LineWidth', 2, 'Color', cols(iBAC, :), 'HandleVisibility', 'off')
    end

    % Plot average FC vs SA:
    scatter(FC, SA, 100, cols(iBAC,:), 'filled',  'Marker', 's')
    
    % Save legend:
    lgd{iBAC} = sprintf('BAC=%u', BAC(iBAC));
end

% Plot settings:
xlabel('Fitness cost', 'FontSize', 17, 'Interpreter', 'Latex')
ylabel('Survival advantage', 'FontSize', 17, 'Interpreter', 'Latex')
xlim([0 1])
ylim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1.0])
yticks([0 0.2 0.4 0.6 0.8 1.0])

legend(lgd, 'Location', 'EastOutside', 'FontSize', 15, 'Interpreter', 'Latex', 'EdgeColor', 'None')
title('Selection plane', 'FontSize', 15, 'Interpreter', 'Latex')

ax = gca;

ax.XColor     = 'k';
ax.YColor     = 'k';
ax.LineWidth  = 2;
ax.TickLength = [0.02 0.02];


% ----------------------------------------------------------------------- %
% (3) Second plot (extinction plane):

% Obtain limits for plot:
try 
    xm = min(CI_SF_W(:, 1));
    xM = 1;
catch
    xm = 0.1*min(SF_W);
    xM = 1;
end

try 
    ym = min(CI_SF_T(:, 1));
    yM = 1;
catch
    ym = 0.1*min(SF_T);
    yM = 1;
end

% ----------------------------------------------------------------------- %
% (3.1) Calculation of the extinction bounds:

% Approximation of the saturation time:
tsat        = log(K/X_0)/((X_W0/X_0)*mu_W + (X_T0/X_0)*mu_T);

% Bound on SF_W for extinction of W in isolation:
ext_W_isol = max(exp(- mu_W*t_g)/D, X_e/(D*K));

% Bound on SF_T for extinction of T in isolation:    
ext_T_isol = max(exp(- mu_T*t_g)/D, X_e/(D*K));

% Bounds for extinction of W in competition:
ext_W      = (exp(- mu_W*tsat)*X_e)/(X_W0*D);

% Bounds for extinction of T in competition:     
ext_T      = (exp(- mu_T*tsat)*X_e)/(X_T0*D);

figure
set(gcf, 'Color', 'w')
hold on
set(gca, 'TickLabelInterpreter', 'Latex', 'XScale', 'log', 'YScale', 'log', 'FontSize', 25)

% ----------------------------------------------------------------------- %
% (3.2) Plot lines delimiting zones:

% Actualise plot limits if the extinction bounds are out of range:
xm = min([xm;0.1*ext_W_isol]);
xM = max([xM;10*ext_W_isol;10*ext_W]);
xM = min(1, xM);
ym = min([ym;0.1*ext_T_isol]);
yM = max([yM;10*ext_T_isol;10*ext_T]);
yM = min(1, yM);

% Auxiliary arrays covering the x- and y-axis:
x = linspace(xm, xM, 100);
y = linspace(ym, yM, 100);

% Plot lines representing the bounds in isolation:
plot(x, ext_T_isol*ones(size(x)), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')
plot(ext_W_isol*ones(size(y)), y, '-k', 'LineWidth', 2, 'HandleVisibility', 'off')

% ----------------------------------------------------------------------- %
% (3.3) Fill regions of the extinction plane:

% (3.3.a) Regions for extinction in isolation:
xaux  = ext_W_isol;
yaux  = ext_T_isol;

% Fill regions:
fill([xm xaux xaux xm], [max(ext_T,yaux) max(ext_T,yaux) yM yM], col_ext_W, 'EdgeColor', 'None',... 
     'FaceAlpha', 0.2, 'LineWidth', 2)   
fill([max(ext_W,xaux) xM xM max(ext_W, xaux)], [ym ym yaux yaux], col_ext_T, 'EdgeColor', 'None',... 
     'FaceAlpha', 0.2, 'LineWidth', 2) 

% Total extinction in isolation:
fill([xm xaux xaux xm], [ym ym yaux yaux], col_ext, 'EdgeColor', 'None',... 
     'FaceAlpha', 0.2, 'LineWidth', 2) 

% (3.3.b) Regions for extinction in competition:
if ext_W > ext_W_isol && ext_T > ext_T_isol                            % Case 1: Region for extinction in competition out of region 
                                                                           % for extinction in isolation;
    % Create legend:
    lgd = {'Ext. W (isol.)', 'Ext. T (isol.)', 'Total Ext. (isol.)', 'Ext. W (comp)', 'Ext. T (comp)', 'Total Ext. (comp.)'};

    xaux = ext_W;
    yaux = ext_T;
    
    % Plot bounds for extinction in competition:
    plot(x, ext_T*ones(size(x)), '--k', 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(ext_W*ones(size(y)), y, '--k', 'LineWidth', 2, 'HandleVisibility', 'off')

    % Extinction in competition:
    fill([ext_W_isol xaux xaux ext_W_isol], [yaux yaux yM yM], col_ext_W, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.6, 'LineWidth', 2)   
    fill([xaux xM xM xaux], [ext_T_isol ext_T_isol yaux yaux], col_ext_T, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.6, 'LineWidth', 2)   

    % Total extinction in competition:
    fill([xm xaux xaux xm], [ext_T_isol ext_T_isol yaux yaux], col_ext, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2) 
    fill([ext_W_isol xaux xaux ext_W_isol], [ym ym ext_T_isol ext_T_isol], col_ext, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2)

elseif ext_W > ext_W_isol                                                % Case 2: Region for extinction in competition of W out of region 
                                                                           % for extinction in isolation of W;
    % Create legend:
    lgd = {'Ext. W (isol.)', 'Ext. T (isol.)', 'Total Ext. (isol.)', 'Ext. W (comp)', 'Total Ext. (comp.)'};
    
    xaux = ext_W;
    yaux = ext_T_isol;
    
    % Plot bound for extinction of W in competition:
    plot(ext_W*ones(size(y)), y, '--k', 'LineWidth', 2, 'HandleVisibility', 'off')

    % Extinction of W in competition:
    fill([ext_W_isol xaux xaux ext_W_isol], [yaux yaux yM yM], col_ext_W, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.6, 'LineWidth', 2)   

    % Total extinction in competition:
    fill([ext_W_isol xaux xaux ext_W_isol], [ym ym yaux yaux], col_ext, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2)

elseif ext_T > ext_T_isol                                                % Case 3: Region for extinction in competition of T out of region 
                                                                           % for extinction in isolation of T;
    % Create legend:
    lgd = {'Ext. W (isol.)', 'Ext. T (isol.)', 'Total Ext. (isol.)', 'Ext. T (comp)', 'Total Ext. (comp.)'};
    xaux = ext_W_isol;
    yaux = ext_T;
    
    % Plot bound for extinction of T in competition:
    plot(x, ext_T*ones(size(x)), '--k', 'LineWidth', 2, 'HandleVisibility', 'off')

    % Extinction in competition:
    fill([xaux xM xM xaux], [ext_T_isol ext_T_isol yaux yaux], col_ext_T, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.6, 'LineWidth', 2)   

    % Total extinction in competition:
    fill([xm xaux xaux xm], [ext_T_isol ext_T_isol yaux yaux], col_ext, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2) 
else
    % Create legend:
    lgd = {'Ext. W (isol.)', 'Ext. T (isol.)', 'Total Ext.(isol.)'};
    
    xaux = ext_W_isol;
    yaux = ext_T_isol;
end

% ----------------------------------------------------------------------- %
% (3.4) Plot selection range:

% intersection point between the equilibrium line and the bound of W:
yint_W = xaux^(1-FC)/(D^FC);

% intersection point between the equilibrium line and the bound of T:
xint_T = (yaux*D^FC)^(1/(1-FC)); 

% Plot line of selective equilibrium:
xint = max(xaux, xint_T);
xx   = linspace(xint, xM, 100);
plot(xx, xx.^(1 - FC)/(D^FC), '-k', 'LineWidth', 2, 'HandleVisibility', 'off')   

% Intersection between the equilibrium line and the maximum x,y-axis:
xint_M  = (yM*D^FC)^(1/(1-FC));
yint_M  = xM^(1-FC)/(D^FC);

% Fill region of selection range (he shape of the range is different, i.e.,
% the number of vertices changes, depending on whether the equilibrium line
% intersects with the bound of T or with the bound of W):
if xint == xaux
    if yint_M < yM
        fill([xaux xM xM xaux], [yaux yaux yint_M yint_W], col_SR_W, 'EdgeColor', 'None',... 
             'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off') 
        fill([xaux xM xM xaux], [yint_W yint_M yM yM], col_SR_T, 'EdgeColor', 'None',... 
             'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    elseif yint_M == yM
        fill([xaux xM xM xaux], [yaux yaux yM yint_W], col_SR_W, 'EdgeColor', 'None',... 
             'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off') 
        fill([xaux xM xaux], [yint_W yM yM], col_SR_T, 'EdgeColor', 'None',... 
             'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    else
        fill([xaux xM xM xint_M xaux], [yaux yaux yM yM yint_W], col_SR_W, 'EdgeColor', 'None',... 
             'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off') 
        fill([xaux xint_M xaux], [yint_W yM yM], col_SR_T, 'EdgeColor', 'None',... 
             'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    end
else  
    if yint_M <= yM 
        fill([xint xM xM], [yaux yaux yint_M], col_SR_W, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off') 
        fill([xaux xint xM xM xaux], [yaux yaux yint_M yM yM], col_SR_T, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    elseif yint_M == yM
        fill([xint xM xM], [yaux yaux yM], col_SR_W, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off') 
        fill([xaux xint xM xaux], [yaux yaux yM yM], col_SR_T, 'EdgeColor', 'None',...
         'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    else
        fill([xint xM xM xint_M], [yaux yaux yM yM], col_SR_W, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off') 
        fill([xaux xint xint_M xaux], [yaux yaux yM yM], col_SR_T, 'EdgeColor', 'None',... 
         'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    end
end

% Plot text indicating selection:
text(0.01*xM, 10*yaux,  'W is selected', 'FontSize', 15, 'Interpreter', 'Latex')
text(5*xaux, 0.5*yM, 'T is selected', 'FontSize', 15, 'Interpreter', 'Latex')
text(1.5*xint, 2*(xint^(1-FC)/(D^FC)), 'Equilibrium', 'FontSize', 15, 'Interpreter','latex', 'Rotation', 30)

% ----------------------------------------------------------------------- %
% (3.5) Plot boxes for each BAC concentration:
for iBAC = 1:nBAC
    scatter(SF_W(iBAC), SF_T(iBAC), 100, cols(iBAC,:), 'filled',  'Marker', 's')

    if numel(CI_SF_W) == numel(CI_SF_T) && numel(CI_SF_W) == 2*nBAC
        fill([CI_SF_W(iBAC, 1) CI_SF_W(iBAC, 2) CI_SF_W(iBAC, 2) CI_SF_W(iBAC, 1)],...
             [CI_SF_T(iBAC, 1) CI_SF_T(iBAC, 1) CI_SF_T(iBAC, 2) CI_SF_T(iBAC, 2)], cols(iBAC, :),...
             'EdgeColor', cols(iBAC, :), 'FaceAlpha', 0.4, 'LineWidth', 2, 'HandleVisibility', 'off')
    end
end

% ----------------------------------------------------------------------- %
% Settings for the figure:

hold off

%Labels for the axes:
xlabel('Survival fraction of W', 'Interpreter', 'Latex', 'FontSize', 17)
ylabel('Survival fraction of T', 'Interpreter','Latex','FontSize', 17);

% Set limits for plot:
xlim([xm xM])
ylim([ym yM])

% Set ticks for axis:
xticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1])
yticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1])

ax = gca;

ax.XAxis.MinorTickValues = [1e-11 1e-9 1e-7 1e-5 1e-3 1e-1];
ax.YAxis.MinorTickValues = [1e-11 1e-9 1e-7 1e-5 1e-3 1e-1];

% Bring axis to front:
set(ax ,'Layer', 'Top')


ax = gca;

ax.XColor     = 'k';
ax.YColor     = 'k';
ax.LineWidth  = 2;
ax.TickLength = [0.02 0.02]; 

% Show legend:
legend(lgd, 'Location', 'NorthOutside', 'Orientation', 'Horizontal', 'FontSize', 12, 'Interpreter', 'Latex', 'EdgeColor', 'None')

%Title for the plot:
title('Extinction plane', 'FontSize', 15, 'Interpreter', 'Latex')

end