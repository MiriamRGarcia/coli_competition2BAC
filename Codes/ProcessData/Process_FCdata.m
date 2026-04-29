%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process_FCdata: Process flow cytometry (FC) data from the competition
% experiment between: 
% between strains: 
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
function Data = Process_FCdata(filename, well_name, ic, plot_res) %#ok<*STOUT>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
% filename = Name of the well to process FC data;
% ic       = Cycle number;
% plot_res = Plot results of data processing (=1) or not (=0);
%
% OUTPUT:
% Data     = Cleaned dataset;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
% (1) Load raw data:

% Call fca_readfcs:
[fcsdat, ~, ~, ~] = fca_readfcs(filename); 

% Size of the dataset:
nData = size(fcsdat, 1);

% Name of the relevant variables:
vars  = {'FSC_H', 'FSC_A', 'SSC_H', 'SSC_A', 'FL1_H', 'FL1_A', 'FL2_H', 'FL2_A',...
         'FL3_H', 'FL3_A', 'FL4_H', 'FL4_A', 'FL5_H', 'FL5_A', 'FL6_H', 'FL6_A', 'FSC_W'};

% Number of relevant variables:
nvars = size(vars, 2);

% Almacenate raw data:
for ivar = 1:nvars
    eval(sprintf('RawData.%s = fcsdat(%u:%u, %u);', vars{ivar}, 1, nData, ivar));
end

% Size of the raw dataset:
nRawData = numel(RawData.FSC_A);

% Find spurious events (debrits, noise) with zero or negative values
% of the variables:
iOut_zero  = find(RawData.FSC_A <= 0 | RawData.SSC_A <= 0 | RawData.FL1_H <= 0 | RawData.FL5_H <= 0);
inOut_zero = (1:nRawData).';
inOut_zero(iOut_zero) = [];

% Actualise data to remove spurious events:
for ivar = 1:nvars
    eval(sprintf('aux = RawData.%s;', vars{ivar}))
    aux = aux(inOut_zero);
    eval(sprintf('Data.%s = aux;', vars{ivar}));
end

% New size of the dataset:
nData = numel(inOut_zero);

% ----------------------------------------------------------------------- %
% (2) First gate to remove noise based on FSC_A vs SSC_A:

% Remove outliers of FSC_A vs SSC_A based on 2d-Mahalanobis distance:
FSC_A     = Data.FSC_A;
SSC_A     = Data.SSC_A;
Z         = [log10(1 + FSC_A) log10(1 + SSC_A)];
d2        = mahal(Z, Z);
threshold = chi2inv(0.99, 2);

% Indexes of non-outliers of FSC_A vs SSC_A:
inOut_FSCvsSSC = find(d2 < threshold);

% Indexes of outliers of FSC_H vs FSC_A:
iOut_FSCvsSSC  = (1:nData).';
iOut_FSCvsSSC(inOut_FSCvsSSC) = [];

% Actualise data to remove outliers of the 2D distribution:
for ivar = 1:nvars
    eval(sprintf('aux = Data.%s;', vars{ivar}))
    aux = aux(inOut_FSCvsSSC);
    eval(sprintf('Data.%s = aux;', vars{ivar}));
end

% New size of the dataset:
nData = numel(inOut_FSCvsSSC);

% ----------------------------------------------------------------------- %
% (2) Remove doublets from SSC-H vs SSC-A:

% Fit line to data:
SSC_H       = Data.SSC_H;
SSC_A       = Data.SSC_A;
p           = polyfit(SSC_H, SSC_A, 1);
SSC_A_fit   = polyval(p, SSC_H);

% Keep 99% of points closest to singlets line:
dist        = abs(SSC_A - SSC_A_fit)./sqrt(1 + p(1)^2);
threshold   = prctile(dist, 99);
inOut_SSC   = find(dist < threshold);

% Indexes of outliers:
iOut_SSC            = (1:nData).';
iOut_SSC(inOut_SSC) = [];

% Actualise data array to exclude doublets:
for ivar = 1:nvars
    eval(sprintf('aux = Data.%s;', vars{ivar}));

    aux = aux(inOut_SSC);
    eval(sprintf('Data.%s = aux;', vars{ivar}));
end

% Actualise size of the dataset:
nData = numel(inOut_SSC);

fprintf('\n >> MESSAGE FROM Process_FCdata for well: %s, cycle: %u', well_name, ic)
fprintf('\n >> %u points were classified as noise, and removed from the dataset.', numel(iOut_FSCvsSSC) + numel(iOut_zero))
fprintf('\n >> %u points were classified as doublets from SSC-H vs SSC-A, and removed from the dataset.', numel(iOut_SSC))
fprintf('\n >> Total percentage of points removed: %.2f%%\n', 100*(nRawData - nData)/nRawData)

% ----------------------------------------------------------------------- %
% (3) Plot results:

if plot_res > 0
    
    % Indexes of all excluded points from raw dataset:
    inOut = ismember(RawData.FSC_H, Data.FSC_H);
    iOut  = (1:nRawData).';
    iOut(inOut) = [];
    
    figure
    set(gcf, 'Color', 'w', 'Position', [674.3333  107.6667  560.0000  420.0000])

    % FSC-H vs FSC-A:
    subplot(2, 3, 1)
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    dscatter(RawData.FSC_A(inOut), RawData.SSC_A(inOut))
    scatter(RawData.FSC_A(iOut), RawData.SSC_A(iOut), 2, 'filled', 'r')
    xlabel('FSC-A', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('SSC-A', 'Interpreter', 'latex', 'FontSize', 12)
    hold off
    title('Raw data', 'Interpreter', 'latex', 'FontSize', 12)
    
    % FSC-H vs FSC-A:
    subplot(2, 3, 4)
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    dscatter(Data.FSC_A, Data.SSC_A)
    xlabel('FSC-A', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('SSC-A', 'Interpreter', 'latex', 'FontSize', 12)
    hold off
    title('Cleaned dataset', 'Interpreter', 'latex', 'FontSize', 12)

    % SSC-H vs SSC-A:
    subplot(2, 3, 2)
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    dscatter(RawData.SSC_H(inOut), RawData.SSC_A(inOut))
    scatter(RawData.SSC_H(iOut), RawData.SSC_A(iOut), 2, 'filled', 'r')
    xlabel('SSC-H', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('SSC-A', 'Interpreter', 'latex', 'FontSize', 12)
    hold off
    title('Raw data', 'Interpreter', 'latex', 'FontSize', 12)
    
    % SSC-H vs SSC-A:
    subplot(2, 3, 5)
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    dscatter(Data.SSC_H, Data.SSC_A)
    xlabel('SSC-H', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('SSC-A', 'Interpreter', 'latex', 'FontSize', 12)
    hold off
    title('Cleaned dataset', 'Interpreter', 'latex', 'FontSize', 12)
    
    % FITC-H vs PE-H:
    subplot(2, 3, 3)
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    dscatter(log10(1 + RawData.FL1_H(inOut)), log10(1 + RawData.FL5_H(inOut)))
    scatter(log10(1 + RawData.FL1_H(iOut)), log10(1 + RawData.FL5_H(iOut)), 2, 'filled', 'r')
    xlabel('FITC-H', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('PE-H', 'Interpreter', 'latex', 'FontSize', 12)
    hold off
    title('Raw data', 'Interpreter', 'latex', 'FontSize', 12)
    
    % FITC-H vs PE-H:
    subplot(2, 3, 6)
    hold on
    set(gca, 'TickLabelInterpreter', 'Latex', 'FontSize', 12)
    dscatter(log10(1 + Data.FL1_H), log10(1 + Data.FL5_H))
    xlabel('FITC-H', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('PE-A', 'Interpreter', 'latex', 'FontSize', 12)
    hold off
    title('Cleaned dataset', 'Interpreter', 'latex', 'FontSize', 12)
    
    pause
    fprintf('\n >> Press any button to continue.')
    close
end

end