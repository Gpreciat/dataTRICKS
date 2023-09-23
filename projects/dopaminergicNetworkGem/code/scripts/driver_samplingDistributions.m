
% Load model and flux distribution
clear
load('~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1/iDopaNeuro1/iDopaNeuro1.mat')
load('~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1/iDopaNeuro1/S_modelRHMC.mat')

pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep ...
    'iDN1', filesep 'iDopaNeuro1'];

% Flux variability
% [minFva, maxFva] = fluxVariability(iDopaNeuro1);

% Reactions to sample
% rxnsToPlot = {'ATPM', 'EX_glc_D[e]'};
% Most informative
ex = find(~cellfun(@isempty, strfind(iDopaNeuro1.rxns, 'EX_')));
dm = find(~cellfun(@isempty, strfind(iDopaNeuro1.rxns, 'DM_')));
sink = find(~cellfun(@isempty, strfind(iDopaNeuro1.rxns, 'sink_')));
external = iDopaNeuro1.rxns([ex; dm; sink]);
externalIdx = find(~cellfun(@isempty, regexp(iDopaNeuro1.rxns([ex; dm; sink]), '(\[e\])')));
external = unique(external(externalIdx));
isIntRxn = ones(length(iDopaNeuro1.rxns), 1);
isIntRxn(findRxnIDs(iDopaNeuro1, external)) = 0;
numOfExchanges = 6;
covS = cov(S_model');
[S_rxns,S_adjusted_norms, S_norms] = getTopKNorms_adaptive(covS, numOfExchanges, isIntRxn);
rxnsToPlot = iDopaNeuro1.rxns(S_rxns);

% Unweighted 2-norm Flux prediction
solution.unWeighted2norm = optimizeCbModel(iDopaNeuro1, 'min', 1e-6);

%Entropy
iDopaNeuro1.osenseStr = 'min';
iDopaNeuro1.cf = 0;
iDopaNeuro1.cr = 0;
iDopaNeuro1.g = 2;
iDopaNeuro1.u0 = 0;
iDopaNeuro1.f = 1;
tcbmParam.method = 'fluxConc';
tcbmParam.printLevel = 0;
tcbmParam.solver = 'mosek';
[solution.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(iDopaNeuro1, tcbmParam);

figure
for i = 1:length(rxnsToPlot)
    subplot(2, 3, i)
    currRxn = findRxnIDs(iDopaNeuro1, rxnsToPlot{i});
    [f, xi] = hist(S_model(currRxn, :), 100);
    f = f / sum(f); % frecuency to probablility
    % Flux distribution
    area(xi, f,'FaceColor', 'r','FaceAlpha', .2, 'EdgeAlpha', .2);
    hold on
    % Mean
    scatter(mean(S_model(currRxn, :)), max(f), 'filled');
    
    % Standard deviation
    errorbar(mean(S_model(currRxn, :)),max(f),std(S_model(currRxn, :)), ...
        'horizontal');

    % minFlux = xline(minFva(currRxn), ':', 'LineWidth', 2);
    % maxFlux = xline(maxFva(currRxn), ':', 'LineWidth', 2);
    % Prediction
    xline(solution.unWeighted2norm.v(currRxn), '-b', 'LineWidth', 2);
    xline(solution.unWeightedTCBMfluxConc.v(currRxn), '--g', 'LineWidth', 2);
    % legend('p1', 'sdD', 'meanD', 'Min flux FVA', 'Max flux FVA', 'prediction', 'location', 'best')
    title(iDopaNeuro1.rxnNames(currRxn), 'FontSize', 18)
    xlabel('Flux distribution', 'FontSize', 16)
    ylabel('Probability', 'FontSize', 16)
    
    % axis
    axis([min([xi solution.unWeighted2norm.x(currRxn)])-((max(xi)-min(xi))*0.1) ...
        max([xi solution.unWeighted2norm.x(currRxn)])+((max(xi)-min(xi))*0.1) ...
        0 max(f)+((max(f)-min(f))*0.1)])
    
    if i == 1
        legend('Marginal distribution', 'Sample Mean', 'Sample standard deviation', '2-norm', 'entropy maximisation', 'location', 'south', 'FontSize', 14)
    end
end

%save the figure
saveas(gcf,[pathSave filesep 'sampleDistributionOfTopPriorityExchangedMetabolites.fig']);
