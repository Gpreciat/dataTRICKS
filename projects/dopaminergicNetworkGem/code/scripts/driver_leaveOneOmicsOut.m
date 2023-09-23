    %% Generate omicsModel model varations for best condition vs best accuracy (Xomics fig)

    clear

    % Define results directory
    pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
        filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

    % Fixed options
    param.printLevel = 1;
    param.setObjective = ''; % No objective function
    feasTol = getCobraSolverParams('LP', 'feasTol');
    param.boundPrecisionLimit = feasTol * 10;
    param.fluxEpsilon = feasTol * 10;
    param.fluxCCmethod = 'fastcc';
    param.curationOverOmics = true;
    param.weightsFromOmics = 1;
    param.metabolomicWeights='mean';
    param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
    param.addCoupledRxns = 1;
    param.nonCoreSinksDemands = 'closeAll';
    param.metabolomicsBeforeModelExtraction = 1;
    param.closeUptakes = true; % Cell culture information
    param.debug = false;

    % Load Recon3D
    inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
        filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
    genericModelName = 'Recon3DModel_301_xomics_input.mat';
    load([inputFolder filesep genericModelName])
    modelGenerationConditions.genericModel.genericModel = model;

    modelConditions = {'All', 'noB', 'noC', 'noM', 'noT'};
    if ~isfile([pathSave filesep 'leaveOneOmicsOutMultidimensionalComparisonStats.mat'])
        % Compare conditions
        models = {'iDopaNeuro1'; 'iDopaNeuro1p1'};
        c = 0;
        for i = 1:length(models)

            % Load model
            load([pathSave filesep models{i} filesep models{i} '.mat'])
            switch models{i}
                case 'iDopaNeuro1'
                    model = iDopaNeuro1;
                case 'iDopaNeuro1p1'
                    model = iDopaNeuro1p1;
            end
            leaveOneOmicsOutModelsDir = [pathSave filesep models{i} filesep 'leaveOneOmicsOutModels'];
            if ~isfolder(leaveOneOmicsOutModelsDir)
                mkdir(leaveOneOmicsOutModelsDir);
            end
            modelGenerationConditions.outputDir = leaveOneOmicsOutModelsDir;

            % specificData
            specificData = model.XomicsToModelSpecificData;
            % Leave one omics out
            modelGenerationConditions.specificData.noB = rmfield(specificData, {'rxns2add', 'rxns2constrain', ...
                'activeReactions', 'activeGenes', 'inactiveGenes'});
            modelGenerationConditions.specificData.noC = rmfield(specificData, 'mediaData');
            modelGenerationConditions.specificData.noM = rmfield(specificData, 'exoMet');
            modelGenerationConditions.specificData.noT = rmfield(specificData, 'transcriptomicData');

            % modelGenerationConditions
            % Note: if the data does not vary it is enough to declare them in param.
            % Here the non-varying conditions are left for demonstrative purposes
            modelGenerationConditions.cobraSolver = {'gurobi'};
            modelGenerationConditions.activeGenesApproach = {model.XomicsToModelParam.activeGenesApproach};
            modelGenerationConditions.transcriptomicThreshold = model.XomicsToModelParam.tresholdT;
            modelGenerationConditions.boundsToRelaxExoMet = {'b'};
            modelGenerationConditions.closeIons = model.XomicsToModelParam.closeIons;
            modelGenerationConditions.tissueSpecificSolver = {model.XomicsToModelParam.tissueSpecificSolver};
            modelGenerationConditions.inactiveGenesTranscriptomics = model.XomicsToModelParam.inactiveGenesTranscriptomics;
            modelGenerationConditions.limitBounds = model.XomicsToModelParam.TolMaxBoundary;

            directories = XomicsToMultipleModels(modelGenerationConditions, param, false);

            % Set parameters
            changeCobraSolver('ibm_cplex','all');
            paramModelTest.tests = 'flux';
            paramModelTest.trainingSet = model.XomicsToModelSpecificData.exoMet;
            paramModelTest.objectives = {'unWeighted0norm'; 'Weighted0normGE'; ...
                    'unWeighted1norm'; 'Weighted1normGE'; 'unWeighted2norm'; ...
                    'Weighted2normGE'; 'unWeightedTCBMflux'; 'unWeightedTCBMfluxConc'};
            paramModelTest.printLevel = 0;

            if i == 1
                % Prepare the table with the accuracy data
                nRows = (1 + length(directories)) * length(models) *  length(paramModelTest.objectives);
                varTypes = {'string', 'string', 'string', 'double', 'double', 'double', 'double', ...
                    'double', 'double', 'double'};
                varNames = {'model', 'condition', 'objective', 'qualitativeBoth', 'spearmanBoth',...
                    'qualitativeModelSec', 'spearmanModelSec', 'qualitativeModelUpt',...
                    'spearmanModelUpt', 'rankOfS'};
                leaveOneOmicsOutMultidimensionalComparisonStats = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
                    'VariableNames', varNames);
            end

            modelConditions = {'All', 'noB', 'noC', 'noM', 'noT'};
            for j = 1:length(modelConditions)

                if j ~= 1
                    load([directories{j - 1} filesep 'Model.mat'])
                    model = Model;
                end

                if isequal(modelConditions{j}, 'noT')
                   model.expressionRxns = zeros(length(model.rxns), 1);
                end

                % test modelPredictiveCapacity
                [comparisonData, summary] = modelPredictiveCapacity(model, paramModelTest);

                for k = 1:length(paramModelTest.objectives)
                    c = c + 1;
                    objectiveBool = ismember(comparisonData.comparisonStats.objective, paramModelTest.objectives{k});
                    bothBool = ismember(comparisonData.comparisonStats.model, "both");
                    modelSecBool = ismember(comparisonData.comparisonStats.model, "modelSec");
                    modelUptBool = ismember(comparisonData.comparisonStats.model, "modelUpt");
                    leaveOneOmicsOutMultidimensionalComparisonStats.model(c) = models{i};
                    leaveOneOmicsOutMultidimensionalComparisonStats.condition(c) = modelConditions{j};
                    leaveOneOmicsOutMultidimensionalComparisonStats.objective(c) = paramModelTest.objectives{k};
                    leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeBoth(c) = comparisonData.comparisonStats.accuracy(objectiveBool & bothBool);
                    leaveOneOmicsOutMultidimensionalComparisonStats.spearmanBoth(c) = comparisonData.comparisonStats.spearman(objectiveBool & bothBool);
                    leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeModelSec(c) = comparisonData.comparisonStats.accuracy(objectiveBool & modelSecBool);
                    leaveOneOmicsOutMultidimensionalComparisonStats.spearmanModelSec(c) = comparisonData.comparisonStats.spearman(objectiveBool & modelSecBool);
                    leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeModelUpt(c) = comparisonData.comparisonStats.accuracy(objectiveBool & modelUptBool);
                    leaveOneOmicsOutMultidimensionalComparisonStats.spearmanModelUpt(c) = comparisonData.comparisonStats.spearman(objectiveBool & modelUptBool);
                    leaveOneOmicsOutMultidimensionalComparisonStats.rankOfS(c) = rank(full(model.S));
                end
            end
        end
        save([pathSave filesep 'leaveOneOmicsOutMultidimensionalComparisonStats.mat'], 'leaveOneOmicsOutMultidimensionalComparisonStats')
    else
        load([pathSave filesep 'leaveOneOmicsOutMultidimensionalComparisonStats.mat'])
    end

    % Set colors and labels
    colors = [0 0.5 1; 0.6392 0.0784 0.1804; 0.4706 0.6706 0.1882; 1 0 1; 1 0.5 0];

    subplot(1, 2, 1)
    iDN1Bool = ismember(leaveOneOmicsOutMultidimensionalComparisonStats.model, "iDopaNeuro1");
    iDN1p1Bool = ismember(leaveOneOmicsOutMultidimensionalComparisonStats.model, "iDopaNeuro1p1");
    for i = 1:length(modelConditions)
        conditionBool = ismember(leaveOneOmicsOutMultidimensionalComparisonStats.condition, modelConditions{i});
        h(i) = text(nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeBoth(conditionBool & iDN1Bool)), ...
            nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.spearmanBoth(conditionBool & iDN1Bool)), ...
            regexprep(['• iDN1\_' modelConditions{i}], '\\_All', ''), 'FontSize', 12,'Color', colors(i, :));
        set(h(i), 'Rotation', 30)
    %     h(i) = text(nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeBoth(conditionBool & iDN1p1Bool)), ...
    %         nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.spearmanBoth(conditionBool & iDN1p1Bool)), ...
    %         ['• iDN1.1\_' modelConditions{i}], 'FontSize', 12,'Color', colors(i, :));
    %     set(h(i), 'Rotation', 30)
    end
    xmin = 0.54;
    xmax = 0.64;
    ymin = 0.30;
    ymax = max(leaveOneOmicsOutMultidimensionalComparisonStats.spearmanBoth,[],'all');
    axis([xmin xmax ymin ymax])
    title('1. Accuracy', 'FontSize', 16)
    xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 12)
    ylabel('Spearman rank', 'fontweight', 'bold', 'FontSize', 12)
    hold on;
    h1 = plot(NaN, NaN, '.', 'Color', colors(1, :));
    h2 = plot(NaN, NaN, '.', 'Color', colors(2, :));
    h3 = plot(NaN, NaN, '.', 'Color', colors(3, :));
    h4 = plot(NaN, NaN, '.', 'Color', colors(4, :));
    h5 = plot(NaN, NaN, '.', 'Color', colors(5, :));
    legend([h1 h2 h3 h4 h5], 'Control model', 'Bibliomics excluded', 'Cell fresh media excluded', 'Exometabolomics excluded', 'Transcriptomics excluded', 'Location', 'best')
    hold off

    subplot(1, 2, 2)
    for i = 1:length(modelConditions)
        conditionBool = ismember(leaveOneOmicsOutMultidimensionalComparisonStats.condition, modelConditions{i});
        h(i) = text(nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeBoth(conditionBool & iDN1Bool)), ...
            nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.rankOfS(conditionBool & iDN1Bool)), ...
            regexprep(['• iDN1\_' modelConditions{i}], '\\_All', ''), 'FontSize', 12,'Color', colors(i, :));
        set(h(i), 'Rotation', 30)
    %     h(i) = text(nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.qualitativeBoth(conditionBool & iDN1p1Bool)), ...
    %         nanmean(leaveOneOmicsOutMultidimensionalComparisonStats.rankOfS(conditionBool & iDN1p1Bool)), ...
    %         ['• iDN1.1\_' modelConditions{i}], 'FontSize', 12,'Color', colors(i, :));
    %     set(h(i), 'Rotation', 30)
    end
    xmin = 0.54;
    xmax = 0.64;
    ymin = min(leaveOneOmicsOutMultidimensionalComparisonStats.rankOfS,[],'all');
    ymax = 1000;
    axis([xmin xmax (ymin - (ymax * 0.05)) ymax])
    title('2. Accuracy vs rank(S)', 'FontSize', 16)
    xlabel('Correct/total predictions', 'fontweight', 'bold', 'FontSize', 12)
    ylabel('rank(S)', 'fontweight', 'bold', 'FontSize', 12)

    savefig([pathSave filesep 'leaveOneOmicsOutMultidimensionalComparison'])
    saveas(gcf,[pathSave filesep 'leaveOneOmicsOutMultidimensionalComparison'],'png')
    saveas(gcf,[pathSave filesep 'leaveOneOmicsOutMultidimensionalComparison'],'eps')