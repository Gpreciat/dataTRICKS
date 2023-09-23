% Identify the most accurate model based on the results of driver_testiDNmodelPredictiveCapacity.m.
% The information gathered from each model is used to identify two sets of conditions:
%   - The conditions for the model with the highest score;
%   - The conditions with the highest score on each group.
%
% The models are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1

clear

format long g
% Define results directory
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

% Load multidimensionalComparisonStats
load([pathSave filesep 'multidimensionalComparisonStats.mat'])
directoriesWithModels = unique(multidimensionalComparisonStats.dirName);

%% Objectives

% Determine the three objectives with the greatest qualitative predictive
% capacity

objectives = unique(multidimensionalComparisonStats.objective);

[qualitativeAccuracy, quantitativeAccuracy] = deal(zeros(size(objectives)));
for i = 1:length(objectives)
    
    % Find rows
    rowsInComparison = strcmp(multidimensionalComparisonStats.objective, ...
        objectives{i});
    
    % Objectives qualitative score
    qualitativeScoreBoth = multidimensionalComparisonStats.qualitativeBoth(rowsInComparison);
    qualitativeScoreModelSec = multidimensionalComparisonStats.qualitativeModelSec(rowsInComparison);
    qualitativeScoreModelUpt = multidimensionalComparisonStats.qualitativeModelUpt(rowsInComparison);
    % Ignores outlier, and NaN values
    toIgnore = isoutlier(qualitativeScoreBoth) | isnan(qualitativeScoreModelSec) | isnan(qualitativeScoreModelUpt);
    qualitativeAccuracy(i) = mean(qualitativeScoreBoth(~toIgnore));
    
    % Objectives quantitative score
    spearmanScoreBoth = multidimensionalComparisonStats.quantitativeBoth(rowsInComparison);
    spearmanScoreModelSec = multidimensionalComparisonStats.quantitativeModelSec(rowsInComparison);
    spearmanScoreModelUpt = multidimensionalComparisonStats.quantitativeModelUpt(rowsInComparison);
    % Ignores outlier, and NaN values
    toIgnore = isoutlier(spearmanScoreBoth) | isnan(spearmanScoreModelSec) | isnan(spearmanScoreModelUpt);
    quantitativeAccuracy(i) = mean(spearmanScoreBoth(~toIgnore));
    
end

% make objectives table
% objectives = objectives;
% table2latex(table(objectives, qualitativeAccuracy, quantitativeAccuracy), 'objectives.tex');

% Save best objectives
qualitativeAccuracy(isnan(qualitativeAccuracy)) = 0;
[~, ia] = sort(qualitativeAccuracy, 'descend');
topObjectives = objectives(ia(1:2));
save([pathSave filesep 'topObjectives'], 'topObjectives')

%% Model

% Determine the model with the greatest qualitative predictive capacity
% based on the obtianed three objectives

[totalQualitativeScoreModels, totalSpearmanScoreModels, noOfRxns, noOfMets, ...
    rankOfS, spearman] = deal(zeros(size(directoriesWithModels)));
for i = 1:length(directoriesWithModels)
    for j = 1:length(topObjectives)
        
        % Find rows
        objectivesBool = ismember(multidimensionalComparisonStats.objective, topObjectives{j});
        dirBool = strcmp(multidimensionalComparisonStats.dirName, directoriesWithModels{i});
        rowsInComparison = dirBool & objectivesBool;
        
        if any(rowsInComparison)
            
            % Qualitative score
            qualitativeScoreBoth = multidimensionalComparisonStats.qualitativeBoth(rowsInComparison);
            qualitativeScoreModelSec = multidimensionalComparisonStats.qualitativeModelSec(rowsInComparison);
            qualitativeScoreModelUpt = multidimensionalComparisonStats.qualitativeModelUpt(rowsInComparison);
            % Ignores outlier, and NaN values
            toIgnore = isoutlier(qualitativeScoreBoth) | isnan(qualitativeScoreModelSec) | isnan(qualitativeScoreModelUpt);
            qualitativeScoreModels(i, j) = mean(qualitativeScoreBoth(~toIgnore));
            
            % Quantitative score
            spearmanScoreBoth = multidimensionalComparisonStats.spearmanBoth(rowsInComparison);
            spearmanScoreModelSec = multidimensionalComparisonStats.spearmanModelSec(rowsInComparison);
            spearmanScoreModelUpt = multidimensionalComparisonStats.spearmanModelUpt(rowsInComparison);
            % Ignores outlier, and NaN values
            toIgnore = isoutlier(spearmanScoreBoth) | isnan(spearmanScoreModelSec) | isnan(spearmanScoreModelUpt);
            spearmanScoreModels(i, j) = mean(spearmanScoreBoth(~toIgnore));
                        
        else
            
            qualitativeScoreModels(i, j) = 0;
            spearmanScoreModels(i, j) = 0;
            
        end
    end
    
    % Number of reactions and metabolites
    noOfMets(i) = unique(multidimensionalComparisonStats.nOfmets(dirBool));
    noOfRxns(i) = unique(multidimensionalComparisonStats.nOfrxns(dirBool));
    rankOfS(i) = unique(multidimensionalComparisonStats.rankOfS(dirBool));
    
    totalQualitativeScoreModels(i) = mean(nonzeros(qualitativeScoreModels(i, :)));
    totalSpearmanScoreModels(i) = mean(nonzeros(spearmanScoreModels(i, :)));
end
% Set zeros to NaN
qualitativeScoreModels(qualitativeScoreModels == 0) = NaN;
spearmanScoreModels(spearmanScoreModels == 0) = NaN;
%% mean difference for each condition

% Separate each condition in the model
directoriesComparison = cell(size(directoriesWithModels, 1), ...
    length(split(directoriesWithModels{1}, '_')));
for i = 1:size(directoriesWithModels, 1)
    directoriesComparison(i, :) = split(directoriesWithModels{i}, '_')';
end

% Identification of data names and groups
c = 0;
fields = cell(size(unique(directoriesComparison), 1), 1);
for i = 1:size(directoriesComparison, 2)
    currentConditions = unique(directoriesComparison(:, i));
    for j = 1:size(currentConditions, 1)
        c = c + 1;
        fields{c, 1} = [currentConditions{j} '_' num2str(i)];
    end
end
dataGroups = regexp(fields, '(\d+)(?!.*\d)', 'match');
dataGroups = [dataGroups{:}]';

% Adjustments to the labels
expression = 'limit|GenesT|Ions|media|Relaxed|transcriptomics';
plotLabels = regexprep(fields, expression, '');
plotLabels = regexprep(plotLabels, 'deleteModelGenes', 'all-rxn');
plotLabels = regexprep(plotLabels, 'oneRxnPerActiveGene', '1-rxn');
plotLabels = regexprep(plotLabels, 'T', '');
plotLabels = regexprep(plotLabels, 'Boundary.100000', '1e5');
plotLabels = regexprep(plotLabels, 'Boundary.10000', '1e4');
plotLabels = regexprep(plotLabels, 'Boundary.1000', '1e3');
plotLabels = regexprep(plotLabels, 'inactive', 'yes');
plotLabels = regexprep(plotLabels, 'NoInactive', 'no');
plotLabels = regexprep(plotLabels, 'Over([^;]*).*?(?=\_)', '');


figure
% Accuracy subplot
subplot(2,3,1);
hold on
limitAx = 0;
for i = 1:length(unique(dataGroups))
    groups{i} = ['group' num2str(i)];
end

% Plot the labels according the mean difference
conditionsAcc = zeros(size(plotLabels));
for i = 1:size(plotLabels, 1)
    parameterLabel = regexprep(fields{i}, '\_(.*)', '');
    groupBool = sum(ismember(directoriesComparison, parameterLabel), 2);
    nanBool = isnan(totalQualitativeScoreModels);
    meanQualitative = (mean(totalQualitativeScoreModels(groupBool & ~nanBool)) - mean(totalQualitativeScoreModels(~nanBool)));
    h1 = text(meanQualitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
        '\_(.*)', '')], 'FontSize', 12,'Color', [0 0.4510 0.7412]);
    set(h1, 'Rotation', 45)
    if limitAx < abs(meanQualitative)
        limitAx = abs(meanQualitative);
    end
    conditionsAcc(i) = meanQualitative;
    if i == size(plotLabels, 1)
        text(0, 0.5, ['  Mean: ' num2str(round(mean(totalQualitativeScoreModels(~isnan(totalQualitativeScoreModels))), 2))], 'FontSize', 14);
    end
end
xline(0, '--r');
axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
hold off
yticks(1:7)
conditionLabels = {'Extraction algorithm', ...
    '#Reactions per active gene', 'Transcript. threshold, log_2(FPKM)', 'Max. flux (uMol/gDW/hr)', ...
    'Transcript. inactive genes', 'Ion exchange', 'Data preference'};
yticklabels(conditionLabels)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12, 'fontweight', 'bold')
title('A. Qualitative accuracy', 'FontSize', 16)
xlabel('Mean difference', 'FontSize', 12, 'fontweight', 'bold')
grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% Euclidean distance subplot
subplot(2,3,2);
hold on
limitAx = 0;
conditionsED = zeros(size(plotLabels));
for i = 1:size(plotLabels, 1)
    % Plot the labels according the mean difference
    parameterLabel = regexprep(fields{i}, '\_(.*)', '');
    groupBool = sum(ismember(directoriesComparison, parameterLabel), 2);
    nanBool = isnan(totalSpearmanScoreModels);
    meanQuantitative = (mean(totalSpearmanScoreModels(groupBool & ~nanBool)) - mean(totalSpearmanScoreModels(~nanBool)));
    h1 = text(meanQuantitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
        '\_(.*)', '')], 'FontSize', 12,'Color', [0.8510 0.3294 0.1020]);
    set(h1, 'Rotation', 45)
    if limitAx < abs(meanQuantitative)
        limitAx = abs(meanQuantitative);
    end
    conditionsED(i) = meanQuantitative;
    if i == size(plotLabels, 1)
        text(0, 0.5, ['  Mean: ' num2str(round(mean(totalSpearmanScoreModels(~nanBool)), 2))...
            ], 'FontSize', 14);
    end
end
xline(0, '--r');
axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
hold off
title('B. Quantitative accuracy', 'FontSize', 16)
xlabel('Mean difference', 'FontSize', 12, 'fontweight', 'bold')
yticks(1:7)
yticklabels({'','', '', '', '', '', ''})
grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% Number of reactions
subplot(2,3,3);
hold on
limitAx = 0;
conditionsRank = zeros(size(plotLabels));
for i = 1:size(plotLabels, 1)
    % Plot the labels according the mean difference
    parameterLabel = regexprep(fields{i}, '\_(.*)', '');
    [idx, ~] = find(ismember(directoriesComparison, parameterLabel));
    meanQuantitative = mean(rankOfS(idx)) - mean(rankOfS);
    h1 = text(meanQuantitative, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
        '\_(.*)', '')], 'FontSize', 12,'Color', [0.9290, 0.6940, 0.1250]);
    set(h1, 'Rotation', 45)
    if limitAx < abs(meanQuantitative)
        limitAx = abs(meanQuantitative);
    end
    conditionsRank(i) = meanQuantitative;
    if i == size(plotLabels, 1)
        text(0, 0.5, ['  Mean: ' num2str(round(mean(rankOfS)))...
            ], 'FontSize', 14);
    end
end
xline(0, '--r');
axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
hold off
title('C. Stoichiometric matrix rank', 'FontSize', 16)
xlabel('Mean difference', 'FontSize', 12, 'fontweight', 'bold')
yticks(1:7)
yticklabels({'', '', '', '', '', '', ''})
grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% savefig([pathSave filesep 'conditions'])
% saveas(gcf,[pathSave filesep 'conditions'],'png')

%% iDopaNeuroCT

% Identify cell-type-specific models
cellTypeSpecificBool = contains(multidimensionalComparisonStats.preferenceCurationOrOmics, 'curationOverOmics');
% Identify the max predictive capacity in cell-type-specific models
maxCTpredictiveCapacity = max(multidimensionalComparisonStats.qualitativeBoth(cellTypeSpecificBool));
% Identify the directory of the most accurate models
iDopaNeuroCTDir = unique(multidimensionalComparisonStats.dirName(multidimensionalComparisonStats.qualitativeBoth == maxCTpredictiveCapacity & cellTypeSpecificBool));
iDopaNeuroCTDir = 'thermoKernel_oneRxnPerActiveGene_transcriptomicsT0_limitBoundary.100000_NoInactiveGenesT_openIons_curationOverOmics';
iDopaNeuroCTBool =  strcmp(multidimensionalComparisonStats.dirName, iDopaNeuroCTDir);
% Select the objectives
objectiveSelected = ismember(multidimensionalComparisonStats.objective, topObjectives{2});
% Predictive capacity
iDopaNeuroCTData(1, 1) = unique(multidimensionalComparisonStats.nOfmets(iDopaNeuroCTBool & objectiveSelected));
iDopaNeuroCTData(2, 1) = unique(multidimensionalComparisonStats.nOfrxns(iDopaNeuroCTBool & objectiveSelected));
iDopaNeuroCTData(3, 1) = unique(multidimensionalComparisonStats.rankOfS(iDopaNeuroCTBool & objectiveSelected));
iDopaNeuroCTData(4, 1) = multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCTBool & objectiveSelected);
iDopaNeuroCTData(5, 1) = multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCTBool & objectiveSelected);
% Coordinates figures E and F
iDNCT_Ecoordinates = [iDopaNeuroCTData(4, 1) iDopaNeuroCTData(5, 1)];
iDNCT_Fcoordinates = [iDopaNeuroCTData(4, 1) iDopaNeuroCTData(3, 1)];

objectivesBool = ismember(multidimensionalComparisonStats.objective, topObjectives);
rowsInComparison = strcmp(multidimensionalComparisonStats.dirName, iDopaNeuroCTDir) & objectivesBool;
iDN1ModelOFValidation = [round(multidimensionalComparisonStats.qualitativeBoth(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelSec(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelUpt(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.spearmanBoth(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelSec(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelUpt(rowsInComparison), 2)'];
% Table
iDopaNeuro1ValidationTable = table(iDN1ModelOFValidation(:, 1), ...
    iDN1ModelOFValidation(:, 2), ...
    ...
    'VariableNames', ...
    multidimensionalComparisonStats.objective(rowsInComparison),...
    'RowNames',...
    {'qualitativeBoth'; ...
    'qualitativeModelSec';...
    'qualitativeModelUpt';...
    'spearmanBoth';...
    'spearmanModelSec';...
    'spearmanModelUpt'});
display(iDopaNeuro1ValidationTable)

% Copy to results folder
modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'multidimensionalModelGeneration'];
if ~isfolder([pathSave filesep 'iDopaNeuroCT'])
    mkdir([pathSave filesep 'iDopaNeuroCT'])
end
copyfile([modelsDir filesep char(iDopaNeuroCTDir)], [pathSave filesep 'iDopaNeuroCT'])
load([pathSave filesep 'iDopaNeuroCT' filesep 'Model'])
iDopaNeuroCT = Model;
save([pathSave filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT'], 'iDopaNeuroCT')
delete([pathSave filesep 'iDopaNeuroCT' filesep 'Model.mat'])
% sbmlModel = writeSBML(iDopaNeuroCT, 'iDopaNeuroCT');

% Print conditions
fid2 = fopen([pathSave filesep 'iDopaNeuroCT' filesep 'conditions.txt'], 'w');
fprintf(fid2, '%s\n', iDopaNeuroCTDir);
fclose(fid2);

%% iDopaNeuroC

% Identify condition-specific models
conditionSpecificBool = contains(multidimensionalComparisonStats.preferenceCurationOrOmics, 'omicsOverCuration');
% Identify the max predictive capacity in condition-specific models
maxCpredictiveCapacity = max(multidimensionalComparisonStats.qualitativeBoth(conditionSpecificBool));
% Identify the directory of the most accurate models
iDopaNeuroCDir = unique(multidimensionalComparisonStats.dirName(multidimensionalComparisonStats.qualitativeBoth == maxCpredictiveCapacity & conditionSpecificBool));
iDopaNeuroCDir = 'thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.100000_NoInactiveGenesT_closedIons_omicsOverCuration';
iDopaNeuroCBool =  strcmp(multidimensionalComparisonStats.dirName, iDopaNeuroCDir);
% Select the objectives
objectiveSelected = ismember(multidimensionalComparisonStats.objective, topObjectives{1});
% Predictive capacity
iDopaNeuroCData(1, 1) = unique(multidimensionalComparisonStats.nOfmets(iDopaNeuroCBool & objectiveSelected));
iDopaNeuroCData(2, 1) = unique(multidimensionalComparisonStats.nOfrxns(iDopaNeuroCBool & objectiveSelected));
iDopaNeuroCData(3, 1) = unique(multidimensionalComparisonStats.rankOfS(iDopaNeuroCBool & objectiveSelected));
iDopaNeuroCData(4, 1) = multidimensionalComparisonStats.qualitativeBoth(iDopaNeuroCBool & objectiveSelected);
iDopaNeuroCData(5, 1) = multidimensionalComparisonStats.spearmanBoth(iDopaNeuroCBool & objectiveSelected);
% Coordinates figures E and F
iDNC_Ecoordinates = [iDopaNeuroCData(4, 1) iDopaNeuroCData(5, 1)];
iDNC_Fcoordinates = [iDopaNeuroCData(4, 1) iDopaNeuroCData(3, 1)];

objectivesBool = ismember(multidimensionalComparisonStats.objective, topObjectives);
rowsInComparison = strcmp(multidimensionalComparisonStats.dirName, iDopaNeuroCDir) & objectivesBool;
iDN1ModelOFValidation = [round(multidimensionalComparisonStats.qualitativeBoth(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelSec(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.qualitativeModelUpt(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.spearmanBoth(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelSec(rowsInComparison), 2)'; ...
    round(multidimensionalComparisonStats.spearmanModelUpt(rowsInComparison), 2)'];
% Table
iDopaNeuro1ValidationTable = table(iDN1ModelOFValidation(:, 1), ...
    iDN1ModelOFValidation(:, 2), ...
    ...
    'VariableNames', ...
    multidimensionalComparisonStats.objective(rowsInComparison),...
    'RowNames',...
    {'qualitativeBoth'; ...
    'qualitativeModelSec';...
    'qualitativeModelUpt';...
    'spearmanBoth';...
    'spearmanModelSec';...
    'spearmanModelUpt'});
display(iDopaNeuro1ValidationTable)

% Copy to results folder
modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'multidimensionalModelGeneration'];
if ~isfolder([pathSave filesep 'iDopaNeuroC'])
    mkdir([pathSave filesep 'iDopaNeuroC'])
end
copyfile([modelsDir filesep char(iDopaNeuroCDir)], [pathSave filesep 'iDopaNeuroC'])
load([pathSave filesep 'iDopaNeuroC' filesep 'Model'])
iDopaNeuroC = Model;
save([pathSave filesep 'iDopaNeuroC' filesep 'iDopaNeuroC'], 'iDopaNeuroC')
delete([pathSave filesep 'iDopaNeuroC' filesep 'Model.mat'])
% sbmlModel = writeSBML(iDopaNeuroC, 'iDopaNeuroC');

% Print conditions
fid2 = fopen([pathSave filesep 'iDopaNeuroC' filesep 'conditions.txt'], 'w');
fprintf(fid2, '%s\n', iDopaNeuroCDir);
fclose(fid2);

%% Models statistics

% Table
bestModelsTable = table(iDopaNeuroCTData, iDopaNeuroCData, ...
    ...
    'VariableNames', ...
    {'iDopaNeuroCT', 'iDopaNeuroC'},...
    'RowNames',...
    {'Metabolites'; ...
    'Reactions';...
    'Rank of S';...
    'Qualitative';...
    'Quantitative (μmol/gDW/hr)'});
display(bestModelsTable)

% Table figure
subplot(2,3,4)
% iDopaNeuro1
text(0.08, 0.6, 'iDopaNeuroCT', 'FontSize', 16);
text(0.18, 0.5, [' ' num2str(iDopaNeuroCTData(1))], 'FontSize', 16);
text(0.18, 0.4, num2str(iDopaNeuroCTData(2)), 'FontSize', 16);
text(0.18, 0.3, [' ' num2str(iDopaNeuroCTData(3))], 'FontSize', 16);
text(0.18, 0.2, [' ' num2str(round(iDopaNeuroCTData(4), 2)) ' '], 'FontSize', 16);
text(0.18, 0.1, num2str(round(iDopaNeuroCTData(5), 2)), 'FontSize', 16);
% iDopaNeuro1.1
text(0.47, 0.6, 'iDopaNeuroC', 'FontSize', 16);
text(0.58, 0.5, num2str(iDopaNeuroCData(1)), 'FontSize', 16);
text(0.58, 0.4, num2str(iDopaNeuroCData(2)), 'FontSize', 16);
text(0.58, 0.3, num2str(iDopaNeuroCData(3)), 'FontSize', 16);
text(0.58, 0.2, [' ' num2str(round(iDopaNeuroCData(4), 2)) ' '], 'FontSize', 16);
text(0.58, 0.1, num2str(round(iDopaNeuroCData(5), 2)), 'FontSize', 16);
yticks(0.1:0.1:0.6)
set(gca,'YTickLabel', {'Quantitative accuracy', 'Qualitative accuracy', 'rank(S)', 'Reactions', 'Metabolites', 'Model'}, 'fontsize', 12, 'fontweight', 'bold')
title('D. Statistics', 'FontSize', 16)
xline(0.4, 'LineWidth',3)
yline(0.55,'LineWidth',3)
axis([0 0.8 0.05 0.65])
set(gca,'XTick',[])


%% 3d Scatter plot

% subplot(2,3,4)
% % consistency, flux, ED
% scatter3(nonzeros(qualitativeScoreModels(:, 1) * 100), log10(quantitativeScoreModels(:, 1)), rankOfS, 'filled')
% hold on
% scatter3(nonzeros(qualitativeScoreModels(:, 2) * 100), log10(quantitativeScoreModels(:, 2)), rankOfS, 'filled', 'r')
% scatter3(nonzeros(qualitativeScoreModels(:, 3) * 100), log10(quantitativeScoreModels(:, 3)), rankOfS, 'filled', 'g')
% hold off
% legend(topObjectives, 'Location', 'best')
% xlabel({'Qualitative: confusion', ' matrix accuracy (%)'}, 'fontweight', 'bold')
% ylabel({'Quantitative:', 'log10(Euclidean norm)'}, 'fontweight', 'bold')
% zlabel({'Rank of the', 'stoichiometric matrix'}, 'fontweight', 'bold')

subplot(2,3,5)
% consistency, flux, ED
scatter(nonzeros(qualitativeScoreModels(:, 1)), spearmanScoreModels(:, 1), 'g')
hold on
scatter(nonzeros(qualitativeScoreModels(:, 2)), spearmanScoreModels(:, 2), 'm')
hold off
title('E. Accuracy', 'FontSize', 16)
xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 12)
ylabel('Spearman rank', 'fontweight', 'bold', 'FontSize', 12)
xmin = 0.55;
xmax = 0.88;
ymin = min(spearmanScoreModels,[],'all');
ymax = max(spearmanScoreModels,[],'all') + max(spearmanScoreModels,[],'all')*0.1;
axis([xmin xmax (ymin - (ymin * 0.1)) (ymax + (ymax * 0.1))])
text(iDNC_Ecoordinates(1), iDNC_Ecoordinates(2), '• iDN\_C', 'FontSize', 14);
text(iDNCT_Ecoordinates(1), iDNCT_Ecoordinates(2), '• iDN\_CT', 'FontSize', 14);
legendLabels = regexprep(topObjectives, 'unWeightedTCBMflux', 'Max. entropy of forward and reverse flux');
legendLabels = regexprep(legendLabels, 'unWeightedTCBMfluxConc', 'Max. entropy of flux and concentration');
legend(legendLabels, 'FontSize', 12)

subplot(2,3,6)
% consistency, flux, ED
scatter(nonzeros(qualitativeScoreModels(:, 1)), rankOfS, 'g');
hold on
scatter(nonzeros(qualitativeScoreModels(:, 2)), rankOfS, 'm')
hold off
title('F. Accuracy vs rank(S)', 'FontSize', 16)
xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 12)
ylabel('rank(S)', 'fontweight', 'bold', 'FontSize', 12)
xmin = 0.55;
xmax = 0.8;
ymin = min(rankOfS,[],'all');
ymax = max(rankOfS,[],'all');
axis([xmin xmax (ymin - (ymin * 0.1)) (ymax + (ymax * 0.1))])
text(iDNC_Fcoordinates(1), iDNC_Fcoordinates(2), '• iDN\_C', 'FontSize', 14);
text(iDNCT_Fcoordinates(1), iDNCT_Fcoordinates(2), '• iDN\_CT', 'FontSize', 14);

% savefig([pathSave filesep 'multidimensionalComparison'])
% saveas(gcf,[pathSave filesep 'multidimensionalComparison'],'png')
% saveas(gcf,[pathSave filesep 'multidimensionalComparison'],'eps')
