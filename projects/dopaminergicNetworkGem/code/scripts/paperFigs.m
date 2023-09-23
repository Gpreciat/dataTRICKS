%% Define paths
clear

% Get a list of all files and folders in this folder.
resultsDir = [pwd filesep];
directoriesWithModels = dir(resultsDir);
directoriesWithModels = struct2cell(directoriesWithModels([directoriesWithModels.isdir]));
directoriesWithModels = directoriesWithModels(1, 3:end)';

%%
[totalFluxAcc, modelSecFluxAcc, modelUptFluxAcc, accuracyRxns, ...
    accuracyMets, totalFluxED, modelSecFluxED, modelUptFluxED] = ...
    deal(zeros(size(directoriesWithModels)));
if 1
    for i = 1:size(directoriesWithModels, 1)
        
        workingDirectory = [resultsDir directoriesWithModels{i} filesep];
        
        if ~isfile([workingDirectory 'accuracy.mat'])
            
            try
                
                load([workingDirectory 'Model.mat'])
                model = Model;
                
                % Identify the core reactions from model.options
                % Active reactions
                activeRxns = [];
                if isfield(model.options, 'setObjective')
                    activeRxns = model.options.setObjective;
                end
                if isfield(model.options, 'rxns2add')
                    activeRxns = [activeRxns; model.options.rxns2add.rxnID];
                end
                if isfield(model.options, 'activeReactions')  && ~isempty(model.options.activeReactions)
                    activeRxns = [activeRxns; model.options.activeReactions];
                end
                if isfield(model.options, 'rxns2constrain') && ~isempty(model.options.rxns2constrain)
                    activeRxns = [activeRxns; model.options.rxns2constrain.rxnID];
                end
                if isfield(model.options, 'coupledRxns') && ~isempty(model.options.coupledRxns)
                    for j = 1:length(model.options.coupledRxns.coupledRxnsList)
                        activeRxns = [activeRxns; split(model.options.coupledRxns.coupledRxnsList{j}, ', ')];
                    end
                end
                if isfield(model.options, 'mediaData') && ~isempty(model.options.mediaData)
                    activeRxns = [activeRxns; model.options.mediaData.rxnID];
                end
                if isfield(model.options, 'exoMet') && ~isempty(model.options.exoMet) && ismember('rxnID', model.options.exoMet.Properties.VariableNames)
                    activeRxns = [activeRxns; model.rxns(ismember(model.rxns, model.options.exoMet.rxnID))];
                end
                activeRxns = unique(activeRxns);
                
                % Inactive reactions
                inactiveRxns = [];
                if isfield(model.options, 'rxns2remove')
                    inactiveRxns = [inactiveRxns; model.options.rxns2remove.rxnID];
                end
                if isfield(model.options, 'inactiveReactions')
                    inactiveRxns = [inactiveRxns; model.options.inactiveReactions];
                end
                inactiveRxns = unique(inactiveRxns);
                
                % Identify the active and inactive genes from model.options
                % Set genes in the correct format
                model.genes = regexprep(model.genes, '\.\d', '');
                model.genes = unique(model.genes);
                model.grRules = regexprep(model.grRules, '\.\d', '');
                if isfield(model,'rules')
                    model = rmfield(model,'rules');
                end
                if isfield(model,'rxnGeneMat')
                    model = rmfield(model,'rxnGeneMat');
                end
                % activeGenes in correct format
                if isfield(model.options, 'activeGenes')
                    if isnumeric(model.options.activeGenes)
                        format long g
                        tmp = cell(length(model.options.activeGenes),1);
                        for j = 1:length(model.options.activeGenes)
                            tmp{j, 1} = num2str(model.options.activeGenes(j));
                        end
                        model.options.activeGenes = tmp;
                    end
                    model.options.activeGenes = regexprep(model.options.activeGenes, '\.\d', '');
                else
                    model.options.activeGenes = [];
                end
                % inactiveGenes in correct format
                if isfield(model.options, 'inactiveGenes')
                    if isnumeric(model.options.inactiveGenes)
                        format long g
                        tmp = cell(length(model.options.inactiveGenes),1);
                        for i = 1:length(model.options.inactiveGenes)
                            tmp{i,1} = num2str(model.options.inactiveGenes(i));
                        end
                        model.options.inactiveGenes = tmp;
                    end
                    model.options.inactiveGenes = regexprep(model.options.inactiveGenes, '\.\d', '');
                else
                    model.options.inactiveGenes = [];
                end
                % transcriptomicData in correct format
                if isfield(model.options, 'transcriptomicData')
                    model.options.transcriptomicData.genes = regexprep(model.options.transcriptomicData.genes, '\.\d', '');
                end
                % proteomicData in correct format
                if isfield(model.options, 'proteomicData')
                    model.options.proteomicData.genes = regexprep(model.options.proteomicData.genes, '\.\d', '');
                end
                
                % Identify active and inactive genes
                inactiveGenes = [];
                if isfield(model.options, 'inactiveGenes')
                    inactiveGenes = [inactiveGenes model.options.inactiveGenes];
                end
                if isfield(model.options, 'transcriptomicData') && ~isempty(model.options.transcriptomicData)
                    if min(model.options.transcriptomicData.expVal) < 0
                        warning('transcriptomic expression in model.options.transcriptomicData.expVal must be non-negative')
                        model.options.transcriptomicData.expVal = model.options.transcriptomicData.expVal + min(model.options.transcriptomicData.expVal);
                    end
                    %identify genes without transcriptomic data using transcript info
                    geneMissingTranscriptomicDataBool = ~ismember(model.genes, model.options.transcriptomicData.genes);
                    % Map transcriptomic genes onto model.genes
                    [bool, locb] = ismember(model.options.transcriptomicData.genes, model.genes);
                    model.geneExpVal(locb(bool)) = model.options.transcriptomicData.expVal(bool);
                    activeModelGeneBool = model.geneExpVal >= exp(model.options.tresholdT);
                    model.options.inactiveGenes = [model.options.inactiveGenes; model.genes(model.geneExpVal < exp(model.options.tresholdT))];
                end
                % Include proteomic data
                if  isfield(model.options, 'proteomicData') && ~isempty(model.options.proteomicData)
                    model.options.proteomicData.Properties.VariableNames = {'genes' 'expVal'};
                    proteomics_data = model.options.proteomicData;
                    proteomics_data.Properties.VariableNames = {'genes' 'expVal'};
                    temp = {};
                    if isnumeric(proteomics_data.genes)
                        for j = 1:length(proteomics_data.genes)
                            temp(end + 1, 1) = {num2str(proteomics_data.genes(j))};
                        end
                        proteomics_data.geneId = temp;
                    else
                        proteomics_data.geneId = proteomics_data.genes;
                    end
                    modelProtein = false(length(proteomics_data.genes),1);
                    for j = 1:length(proteomics_data.geneId)
                        if ismember(proteomics_data.geneId(j), model.genes)
                            modelProtein(j) = 1;
                        end
                    end
                    modelProtein_Entrez = proteomics_data.genes(modelProtein);
                    modelProtein_Value = proteomics_data.expVal(modelProtein);
                    [activeProteins, inactiveProteins] = activeProteinList(modelProtein_Value, modelProtein_Entrez, model.options.tresholdP, model.options.printLevel-1);
                    if isnumeric(activeProteins)
                        for j = 1:length(activeProteins)
                            temp(j) = {num2str(activeProteins(j))};
                        end
                        activeProteins = temp;
                    end
                    for j=1:length(model.genes)
                        if ismember(model.genes(j), activeProteins)
                            activeModelGeneBool(j) = 1;
                        end
                    end
                end
                % Active genes from transcriptomic and proteomic data
                if ~any(activeModelGeneBool)
                    activeEntrezGeneID = [];
                else
                    activeEntrezGeneID = model.genes(activeModelGeneBool);
                end
                % Active genes from manual curation
                if isfield(model.options, 'activeGenes')
                    activeEntrezGeneID = [activeEntrezGeneID; model.options.activeGenes];
                end
                %unique genes
                activeEntrezGeneID = unique(activeEntrezGeneID);
                
                % Add gene data to the list of reactions
                [~, ~, rxnInGenes, ~] = deleteModelGenes(model, activeEntrezGeneID(ismember(activeEntrezGeneID, model.genes)));
                activeRxns = unique([activeRxns; rxnInGenes]);
                [~, ~, rxnInGenes, ~] = deleteModelGenes(model, activeEntrezGeneID(ismember(activeEntrezGeneID, model.genes)));
                inactiveGenes = unique([inactiveGenes; rxnInGenes]);
                
                % Identify the core metabolites from model.options
                % Active metabolites
                activeMets = [];
                if isfield(model.options, 'presentMetabolites')
                    activeMets = model.options.presentMetabolites.mets;
                end
                externalRxnsBool = ~cellfun(@isempty, regexp(activeRxns, 'EX\_|DM\_|sink\_'));
                if any(externalRxnsBool)
                    activeMets = unique([activeMets; unique(regexprep(activeRxns(externalRxnsBool), 'EX\_|DM\_|sink\_', ''))]);
                end
                
                inactiveMets = [];
                if isfield(model.options, 'absentMetabolites')
                    inactiveMets = model.options.absentMetabolites.mets;
                end
                
                % activeInactiveRxn
                totalRxns = [activeRxns; inactiveRxns];
                model.activeInactiveRxn = double(ismember(model.rxns, totalRxns));
                model.activeInactiveRxn(ismember(model.rxns, inactiveRxns)) = -1;
                
                % presentAbsentMet
                totalMets = [activeMets; inactiveMets];
                model.presentAbsentMet = double(ismember(model.mets, totalMets));
                if ~isempty(inactiveMets)
                    model.presentAbsentMet(ismember(model.mets, inactiveMets)) = -1;
                end
                
                [accuracySummary, fullReport, comparisonStats] = ...
                    testModelAccuracy(model, model.presentAbsentMet, ...
                    model.activeInactiveRxn, model.options.exoMet, ...
                    'fluxConsistent', 0);
                save([workingDirectory 'accuracy.mat'], 'accuracySummary', ...
                    'fullReport', 'comparisonStats')
                
                totalFluxAcc(i) = accuracySummary.Value(1);
                modelSecFluxAcc(i) = accuracySummary.Value(2);
                modelUptFluxAcc(i) = accuracySummary.Value(3);
                accuracyRxns(i) = accuracySummary.Value(4);
                accuracyMets(i) = accuracySummary.Value(5);
                totalFluxED(i) = accuracySummary.Value(6);
                modelSecFluxED(i) = accuracySummary.Value(7);
                modelUptFluxED(i) = accuracySummary.Value(8);
                
            catch
            end
        else
            load([workingDirectory 'accuracy.mat'])
            
            totalFluxAcc(i) = accuracySummary.Value(1);
            modelSecFluxAcc(i) = accuracySummary.Value(2);
            modelUptFluxAcc(i) = accuracySummary.Value(3);
            accuracyRxns(i) = accuracySummary.Value(4);
            accuracyMets(i) = accuracySummary.Value(5);
            totalFluxED(i) = accuracySummary.Value(6);
            modelSecFluxED(i) = accuracySummary.Value(7);
            modelUptFluxED(i) = accuracySummary.Value(8);
        end
    end
else
    for i = 1:size(directoriesWithModels, 1)
        
        workingDirectory = [resultsDir directoriesWithModels{i} filesep];
        
        if ~isfile([workingDirectory 'accuracy.mat'])
            
            [accuracySummary, ~, ~] = testModelAccuracy(model, ...
                model.presentAbsentMet, model.activeInactiveRxn, ...
                model.options.exoMet, 'fluxConsistent', 0);
            
            save([workingDirectory 'accuracy.mat'], 'accuracySummary', ...
                'fullReport', 'comparisonStats')
            
            totalFluxAcc(i) = accuracySummary.Value(1);
            modelSecFluxAcc(i) = accuracySummary.Value(2);
            modelUptFluxAcc(i) = accuracySummary.Value(3);
            accuracyRxns(i) = accuracySummary.Value(4);
            accuracyMets(i) = accuracySummary.Value(5);
            totalFluxED(i) = accuracySummary.Value(6);
            modelSecFluxED(i) = accuracySummary.Value(7);
            modelUptFluxED(i) = accuracySummary.Value(8);
            
        else
            load([workingDirectory 'accuracy.mat'])
            
            totalFluxAcc(i) = accuracySummary.Value(1);
            modelSecFluxAcc(i) = accuracySummary.Value(2);
            modelUptFluxAcc(i) = accuracySummary.Value(3);
            accuracyRxns(i) = accuracySummary.Value(4);
            accuracyMets(i) = accuracySummary.Value(5);
            totalFluxED(i) = accuracySummary.Value(6);
            modelSecFluxED(i) = accuracySummary.Value(7);
            modelUptFluxED(i) = accuracySummary.Value(8);
            
        end
    end
end

%% Mean accuracy and euclidean distance

accuracy = mean([modelSecFluxAcc, modelUptFluxAcc, accuracyRxns, accuracyMets], 2);
ED = mean([modelSecFluxED, modelUptFluxED], 2);

%% mean difference for each condition

% Separate each condition in the model
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
plotLabels = regexprep(plotLabels, 'deleteModelGenes', 'dMG');
plotLabels = regexprep(plotLabels, 'oneRxnPerActiveGene', 'oRPAG');
plotLabels = regexprep(plotLabels, 'NoInactive', 'no');
plotLabels = regexprep(plotLabels, 'inactive', 'yes');

figure
% Accuracy subplot
subplot(1,2,1);
hold on
limitAx = 0;
groups = {'solver','inputData', 'activeGenesApproach', 'transcriptomicThreshold',...
    'toRelax', 'ions', 'tissueSpecificSolver', 'BPLexoMet', 'inactiveGenesTranscriptomics',...
    'limitBounds'};
conditionsAcc = zeros(size(plotLabels));
for i = 1:size(plotLabels, 1)
    % Plot the labels according the mean difference
    groupLabel = regexprep(fields{i}, '\_(.*)', '');
    [idx, ~] = find(ismember(directoriesComparison, groupLabel));
    meanDiffAccuracy = mean(accuracy(idx)) - mean(accuracy);
    h1 = text(meanDiffAccuracy, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
        '\_(.*)', '')], 'FontSize', 12,'Color', [0 0.4510 0.7412]);
    set(h1, 'Rotation', 45)
    if limitAx < abs(meanDiffAccuracy)
        limitAx = abs(meanDiffAccuracy);
    end
    conditionsAcc(i) = meanDiffAccuracy;
end
xline(0, '--r');
axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
hold off
title('Accuracy')
xlabel('Mean difference', 'FontSize', 20)
yticks(1:10)
yticklabels({'','', '', '', '', '', '', '', '', ''})
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',20)
grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

% Euclidean distance subplot
subplot(1,2,2);
hold on
limitAx = 0;
conditionsED = zeros(size(plotLabels));
for i = 1:size(plotLabels, 1)
    % Plot the labels according the mean difference
    groupLabel = regexprep(fields{i}, '\_(.*)', '');
    [idx, ~] = find(ismember(directoriesComparison, groupLabel));
    meanDiffED = mean(ED(idx)) - mean(ED);
    h1 = text(meanDiffED, str2double(dataGroups{i}), ['• ' regexprep(plotLabels{i},...
        '\_(.*)', '')], 'FontSize', 12,'Color', [0.8510 0.3294 0.1020]);
    set(h1, 'Rotation', 45)
    if limitAx < abs(meanDiffED)
        limitAx = abs(meanDiffED);
    end
    conditionsED(i) = meanDiffED;
end
xline(0, '--r');
axis([-limitAx limitAx 0 size(unique(dataGroups), 1) + 1])
hold off
title('Euclidean distance')
xlabel('Mean difference', 'FontSize', 20)
yticks(1:10)
yticklabels(groups)
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',20)
% set(gca,'XTickLabelMode','auto')
grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';

%% Best conditions vs best accuracy

bestConditionsDir = [];
uniqueDataGroups = unique(dataGroups);
for i = 1:size(uniqueDataGroups, 1)
    bestConditionsDir = strcat(bestConditionsDir, regexprep(fields(conditionsAcc == max(conditionsAcc(~cellfun(@isempty, regexp(fields, ['\_' num2str(i) '$']))))), '\_(.*)', ''));
    if i ~= size(uniqueDataGroups, 1)
        bestConditionsDir = strcat(bestConditionsDir, '_');
    end
end
bestAccuracyDir = directoriesWithModels(accuracy == max(accuracy));

%% Generate omicsModel model varations for best condition vs best accuracy (Xomics fig)

load([resultsDir 'initialDataXomics_inputData.mat']);

if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end

% Data type
options.multipleModels.dataType = {'noB', 'noC', 'noM', 'noT'};

for i = 1:2
    
    % directory to save the models
    if i == 1
        options.multipleModels.outputDir = strcat(resultsDir, bestAccuracyDir);
        conditions = strsplit(char(bestAccuracyDir), '\_');
    else
        options.multipleModels.outputDir = strcat(resultsDir, bestConditionsDir);
        conditions = strsplit(char(bestConditionsDir), '\_');
    end
    
    % Cobra solver
    options.multipleModels.cobraSolver = conditions(1);
    % Input data
    options.multipleModels.inputData = conditions(2);
    % Active genes approach
    options.multipleModels.activeGenesApproach = conditions(3);
    % Transcriptomic threshold
    options.multipleModels.transcriptomicThreshold = conditions(4);
    % Bounds to relax
    options.multipleModels.toRelax = conditions(5);
    % Ions exchange
    options.multipleModels.ions = conditions(6);
    % Tissue specific solver
    options.multipleModels.tissueSpecificSolver = conditions(7);
    % Force activation
    options.multipleModels.forceActivation = conditions(8);
    % Use inactive genes from transcriptomics
    options.multipleModels.inactiveGenesTranscriptomics = conditions(9);
    % Limit bounds
    options.multipleModels.limitBounds = conditions(10);
    
    XomicsToMultipleModels(model, options);
    
end

%% test accuracy

figure
for i = 1:2
    
    % directory to save the models
    if i == 1
        dirToCheck = strcat(resultsDir, bestAccuracyDir);
    else
        dirToCheck = strcat(resultsDir, bestConditionsDir);
    end
    
    models = dir([char(dirToCheck) filesep '*Model.mat']);
    models = {models.name};
    
    for j = 1:length(models)
        
        eval('load([char(dirToCheck) filesep models{j}])')
        eval(sprintf('model = %s;', regexprep(models{j}, '.mat', '')))
        
        if ~isfield(model, 'activeInactiveRxn')
            
            % Identify the core reactions from options
            % Active reactions
            activeRxns = [];
            if isfield(options, 'setObjective')
                activeRxns = options.setObjective;
            end
            if isfield(options, 'rxns2add')
                activeRxns = [activeRxns; options.rxns2add.rxnID];
            end
            if isfield(options, 'activeReactions')  && ~isempty(options.activeReactions)
                activeRxns = [activeRxns; options.activeReactions];
            end
            if isfield(options, 'rxns2constrain') && ~isempty(options.rxns2constrain)
                activeRxns = [activeRxns; options.rxns2constrain.rxnID];
            end
            if isfield(options, 'coupledRxns') && ~isempty(options.coupledRxns)
                for k = 1:length(options.coupledRxns.coupledRxnsList)
                    activeRxns = [activeRxns; split(options.coupledRxns.coupledRxnsList{k}, ', ')];
                end
            end
            if isfield(options, 'mediaData') && ~isempty(options.mediaData)
                activeRxns = [activeRxns; options.mediaData.rxnID];
            end
            if isfield(options, 'exoMet') && ~isempty(options.exoMet) && ismember('rxnID', options.exoMet.Properties.VariableNames)
                activeRxns = [activeRxns; model.rxns(ismember(model.rxns, options.exoMet.rxnID))];
            end
            activeRxns = unique(activeRxns);
            
            % Inactive reactions
            inactiveRxns = [];
            if isfield(options, 'rxns2remove')
                inactiveRxns = [inactiveRxns; options.rxns2remove.rxnID];
            end
            if isfield(options, 'inactiveReactions')
                inactiveRxns = [inactiveRxns; options.inactiveReactions];
            end
            inactiveRxns = unique(inactiveRxns);
            
            % Identify the active and inactive genes from options
            % Set genes in the correct format
            model.genes = regexprep(model.genes, '\.\d', '');
            model.genes = unique(model.genes);
            model.grRules = regexprep(model.grRules, '\.\d', '');
            if isfield(model,'rules')
                model = rmfield(model,'rules');
            end
            if isfield(model,'rxnGeneMat')
                model = rmfield(model,'rxnGeneMat');
            end
            % activeGenes in correct format
            if isfield(options, 'activeGenes')
                if isnumeric(options.activeGenes)
                    format long g
                    tmp = cell(length(options.activeGenes),1);
                    for k = 1:length(options.activeGenes)
                        tmp{k, 1} = num2str(options.activeGenes(k));
                    end
                    options.activeGenes = tmp;
                end
                options.activeGenes = regexprep(options.activeGenes, '\.\d', '');
            else
                options.activeGenes = [];
            end
            % inactiveGenes in correct format
            if isfield(options, 'inactiveGenes')
                if isnumeric(options.inactiveGenes)
                    format long g
                    tmp = cell(length(options.inactiveGenes),1);
                    for k = 1:length(options.inactiveGenes)
                        tmp{k,1} = num2str(options.inactiveGenes(k));
                    end
                    options.inactiveGenes = tmp;
                end
                options.inactiveGenes = regexprep(options.inactiveGenes, '\.\d', '');
            else
                options.inactiveGenes = [];
            end
            % transcriptomicData in correct format
            if isfield(options, 'transcriptomicData')
                options.transcriptomicData.genes = regexprep(options.transcriptomicData.genes, '\.\d', '');
            end
            % proteomicData in correct format
            if isfield(options, 'proteomicData')
                options.proteomicData.genes = regexprep(options.proteomicData.genes, '\.\d', '');
            end
            
            % Identify active and inactive genes
            inactiveGenes = [];
            if isfield(options, 'inactiveGenes')
                inactiveGenes = [inactiveGenes options.inactiveGenes];
            end
            if isfield(options, 'transcriptomicData') && ~isempty(options.transcriptomicData)
                if min(options.transcriptomicData.expVal) < 0
                    warning('transcriptomic expression in options.transcriptomicData.expVal must be non-negative')
                    options.transcriptomicData.expVal = options.transcriptomicData.expVal + min(options.transcriptomicData.expVal);
                end
                %identify genes without transcriptomic data using transcript info
                geneMissingTranscriptomicDataBool = ~ismember(model.genes, options.transcriptomicData.genes);
                % Map transcriptomic genes onto model.genes
                [bool, locb] = ismember(options.transcriptomicData.genes, model.genes);
                model.geneExpVal(locb(bool)) = options.transcriptomicData.expVal(bool);
                activeModelGeneBool = model.geneExpVal >= exp(model.options.tresholdT);
                options.inactiveGenes = [options.inactiveGenes; model.genes(model.geneExpVal < exp(model.options.tresholdT))];
            end
            % Include proteomic data
            if  isfield(options, 'proteomicData') && ~isempty(options.proteomicData)
                options.proteomicData.Properties.VariableNames = {'genes' 'expVal'};
                proteomics_data = options.proteomicData;
                proteomics_data.Properties.VariableNames = {'genes' 'expVal'};
                temp = {};
                if isnumeric(proteomics_data.genes)
                    for k = 1:length(proteomics_data.genes)
                        temp(end + 1, 1) = {num2str(proteomics_data.genes(k))};
                    end
                    proteomics_data.geneId = temp;
                else
                    proteomics_data.geneId = proteomics_data.genes;
                end
                modelProtein = false(length(proteomics_data.genes),1);
                for k = 1:length(proteomics_data.geneId)
                    if ismember(proteomics_data.geneId(k), model.genes)
                        modelProtein(k) = 1;
                    end
                end
                modelProtein_Entrez = proteomics_data.genes(modelProtein);
                modelProtein_Value = proteomics_data.expVal(modelProtein);
                [activeProteins, inactiveProteins] = activeProteinList(modelProtein_Value, modelProtein_Entrez, options.tresholdP, options.printLevel-1);
                if isnumeric(activeProteins)
                    for k = 1:length(activeProteins)
                        temp(k) = {num2str(activeProteins(k))};
                    end
                    activeProteins = temp;
                end
                for k=1:length(model.genes)
                    if ismember(model.genes(k), activeProteins)
                        activeModelGeneBool(k) = 1;
                    end
                end
            end
            % Active genes from transcriptomic and proteomic data
            if ~any(activeModelGeneBool)
                activeEntrezGeneID = [];
            else
                activeEntrezGeneID = model.genes(activeModelGeneBool);
            end
            % Active genes from manual curation
            if isfield(options, 'activeGenes')
                activeEntrezGeneID = [activeEntrezGeneID; options.activeGenes];
            end
            %unique genes
            activeEntrezGeneID = unique(activeEntrezGeneID);
            
            % Add gene data to the list of reactions
            [~, ~, rxnInGenes, ~] = deleteModelGenes(model, activeEntrezGeneID(ismember(activeEntrezGeneID, model.genes)));
            activeRxns = unique([activeRxns; rxnInGenes]);
            [~, ~, rxnInGenes, ~] = deleteModelGenes(model, activeEntrezGeneID(ismember(activeEntrezGeneID, model.genes)));
            inactiveGenes = unique([inactiveGenes; rxnInGenes]);
            
            % Identify the core metabolites from options
            % Active metabolites
            activeMets = [];
            if isfield(options, 'presentMetabolites')
                activeMets = options.presentMetabolites.mets;
            end
            externalRxnsBool = ~cellfun(@isempty, regexp(activeRxns, 'EX\_|DM\_|sink\_'));
            if any(externalRxnsBool)
                activeMets = unique([activeMets; unique(regexprep(activeRxns(externalRxnsBool), 'EX\_|DM\_|sink\_', ''))]);
            end
            
            inactiveMets = [];
            if isfield(options, 'absentMetabolites')
                inactiveMets = options.absentMetabolites.mets;
            end
            
            % activeInactiveRxn
            totalRxns = [activeRxns; inactiveRxns];
            model.activeInactiveRxn = double(ismember(model.rxns, totalRxns));
            model.activeInactiveRxn(ismember(model.rxns, inactiveRxns)) = -1;
            
            % presentAbsentMet
            totalMets = [activeMets; inactiveMets];
            model.presentAbsentMet = double(ismember(model.mets, totalMets));
            if ~isempty(inactiveMets)
                model.presentAbsentMet(ismember(model.mets, inactiveMets)) = -1;
            end
            
        end
        
        % Test the accuracy for each model
        [accuracySummary, ~, ~] = testModelAccuracy(model, ...
            model.presentAbsentMet, model.activeInactiveRxn, ...
            options.exoMet, 'fluxConsistent', 0);
        comarisonMatrix(j, 1) = accuracySummary.Value(2);
        comarisonMatrix(j, 2) = accuracySummary.Value(3);
        comarisonMatrix(j, 3) = accuracySummary.Value(4);
        comarisonMatrix(j, 4) = accuracySummary.Value(5);
        comarisonMatrix(j, 5) = mean(comarisonMatrix(j, 1:4));
        comarisonMatrix(j, 6) = accuracySummary.Value(7);
        comarisonMatrix(j, 7) = accuracySummary.Value(8);
        comarisonMatrix(j, 8) = size(model.mets, 1);
        comarisonMatrix(j, 9) = size(model.rxns, 1);
        
    end
    
    % Plot comparison    
    subplot(1, 2, i)
    bar(comarisonMatrix(:, 1:4) * 100,'FaceAlpha', .3, 'EdgeAlpha', .5)
    if i == 1
        title('Higher accuracy','FontSize',16)
    else
        title('Comparison','FontSize',16)
    end
    ylabel('Accuracy (%)','FontSize',16)
    legend({'Upt prediction', 'Sec prediction', 'Rxns consistency', ...
        'Mets consistency'}, 'Location', 'northeastoutside')
    xticks(1:length(models))
    xticklabels(regexprep(models, '.mat', ''))
    xtickangle(45)
    
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',regexprep(models, '.mat', ''),'fontsize', 16)
    
    ax = gca;
    text(1:length(models), ones(1, length(models)) * min(ax.YLim) + (max(ax.YLim) * 0.1),...
        strcat('ED=', num2str(round(mean(comarisonMatrix(:, 6:7)'))')),...
        'vert','bottom','horiz','center','FontSize',14, 'fontweight', 'bold');
    text(1:length(models), ones(1, length(models)) * min(ax.YLim) + (max(ax.YLim) * 0.05),...
        strcat('mets:', num2str(comarisonMatrix(:, 8))),...
        'vert','bottom','horiz','center','FontSize',14, 'fontweight', 'bold');
    text(1:length(models), ones(1, length(models)) * min(ax.YLim) + (max(ax.YLim) * 0.0),...
        strcat('rxns:', num2str(comarisonMatrix(:, 9))),...
        'vert','bottom','horiz','center','FontSize',14, 'fontweight', 'bold');
end




