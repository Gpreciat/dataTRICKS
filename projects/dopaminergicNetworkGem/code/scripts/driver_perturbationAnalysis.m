

% The iDopaNeuro model is perturbed by the inhibition of mitochondrial
% complexes I and V. In a cell culture with dopaminergic neurons, the
% corresponding perturbations are made by treating them with rotenone or
% oligomycin to mimic the inhibition of mitochondrial complexes I and V,
% respectively. Finally, the in silico predictions are validated against
% the experimental data.
%
% The results are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1

clear

% Results directory
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

% Load the top objectives
load([pathSave filesep 'topObjectives'])

%% Prepare experimental data

exoMetabolomicsDataDir = ['~' filesep 'work' filesep 'sbgCloud' filesep ...
    'programReconstruction' filesep 'projects' filesep 'exoMetDN' filesep ...
    'data' filesep 'omics'];

% Read the files
exometData = readtable([exoMetabolomicsDataDir filesep ...
    '20210118_validationData_K7_final_ordered_umol_gDW_h.csv']);
perturbationData = readtable([exoMetabolomicsDataDir filesep ...
    '20210118_validationData_K7_finalAll_ordered_umol_gDW_h.csv']);
rxns = cell(size(perturbationData.variable));
[mlt, nlt] = size(exometData.variable);

% Complete table
for i = 1:mlt
    rxns(ismember(perturbationData.variable, exometData.variable{i})) = ...
        exometData.exRxns(i);
end
perturbationData = [perturbationData rxns];
perturbationData.Properties.VariableNames{'Var10'} = 'rxnID';
perturbationData.Properties.VariableNames{'sds'} = 'SD';
perturbationData.Properties.VariableNames{'variable'} = 'name';
perturbationData.Properties.VariableNames{'platform'} = 'Platform';

% c1TrainingSet
rotenoneBool = strcmp(perturbationData.condition, 'rotenone');
C1TrainingSet = perturbationData(rotenoneBool, :);
C1TrainingSet = sortrows(C1TrainingSet, 'rxnID');
% c5TrainingSet
rotenoneBool = strcmp(perturbationData.condition, 'oligomycin');
C5TrainingSet = perturbationData(rotenoneBool, :);
C5TrainingSet = sortrows(C5TrainingSet, 'rxnID');

%% In silico vs in vitro

models = {'iDopaNeuro1Core'; 'iDopaNeuro1ConditionType'};
for j = 1:length(models)
    
    % Load model
    load([pathSave filesep models{j} filesep models{j} '.mat'])
    switch models{j}
            case 'iDopaNeuro1Core'
                ctrl = iDopaNeuro1Core;
            case 'iDopaNeuro1ConditionType'
                ctrl = iDopaNeuro1ConditionType;
    end
    
    % In silico perturbations
    % Mitochondrial Complex I inhibition 100%
    C1 = changeRxnBounds(ctrl, 'NADH2_u10m', 0, 'b');
    % Mitochondrial Complex V inhibition 100%
    C5 = changeRxnBounds(ctrl, 'ATPS4m', 0, 'b');
    % GBA deletion 100%
    gbaIdx = ctrl.genes(strmatch('2629.', ctrl.genes));
    [GBA, ~, ~, ~] = deleteModelGenes(ctrl, gbaIdx, 0);
        
    %% Validate predictions with the best objective functions
    
    % ctrlTrainingSet
    ctrlTrainingSet = ctrl.XomicsToModelSpecificData.exoMet;
    rowsToDeleteBool = ismember(ctrlTrainingSet.rxnID, C1TrainingSet.rxnID);
    ctrlTrainingSet = ctrlTrainingSet(rowsToDeleteBool, :);
    ctrlTrainingSet = sortrows(ctrlTrainingSet, 'rxnID');
    
    % Compare predictions
    perturbedModels = {'ctrl'; 'C1'; 'C5'; 'GBA'};
    for i = 1:length(perturbedModels)
        
        % Set parameters
        param.tests = 'flux';
        param.objectives = topObjectives(2);
        param.printLevel = 0;
        % Select model and trainingSet
        switch perturbedModels{i}
            case 'ctrl'
                modelTest = ctrl;
                trainingSet = ctrlTrainingSet;
                trainingSet(~ismember(trainingSet.rxnID, modelTest.rxns), :) = [];
                ctrlTrainingSet = trainingSet;
            case 'C1'
                modelTest = C1;
                trainingSet = C1TrainingSet;
                trainingSet(~ismember(trainingSet.rxnID, modelTest.rxns), :) = [];
                C1TrainingSet = trainingSet;
            case 'C5'
                modelTest = C5;
                trainingSet = C5TrainingSet;
                trainingSet(~ismember(trainingSet.rxnID, modelTest.rxns), :) = [];
                C5TrainingSet = trainingSet;
            case 'GBA'
                solutions = modelMultipleObjectives(GBA, param);
                continue
        end
        trainingSet(~ismember(trainingSet.rxnID, modelTest.rxns), :) = [];
        param.trainingSet = trainingSet;
        
        % test modelPredictiveCapacity
        [comparisonData, summary] = modelPredictiveCapacity(modelTest, param);
        comparisonData.fullReport.data = repmat("Model", size(comparisonData.fullReport.model));
        perturbedModelsComparison.(perturbedModels{i}) = comparisonData;
        
    end
    
    % Predicted values
    plotParam.fullReport.ctrl = perturbedModelsComparison.ctrl.fullReport;
    plotParam.fullReport.C1 = perturbedModelsComparison.C1.fullReport;
    plotParam.fullReport.C5 = perturbedModelsComparison.C5.fullReport;
    
    % Experimental values
    plotParam.vExp.ctrl = ctrlTrainingSet;
    plotParam.vExp.C1 = C1TrainingSet;
    plotParam.vExp.C5 = C5TrainingSet;
    
    % Plot data
    plotParam.data = "Model";
    plotParam.objective = topObjectives{2}; 
    plotParam.labelType = 'metabolite';
    plotParam.saveFigures = 1;
    plotParam.savePath = [pathSave filesep models{j}];
    
    plotExperimentalvsPredictedUptakeSecretion(plotParam)
    display(models{j})
    
    % table
    for i = 1:length(perturbedModels(1:3))
       columnsTable(:, i) = [round(100 * perturbedModelsComparison.(perturbedModels{i}).comparisonStats.accuracy(1)), ...
           round(100 * perturbedModelsComparison.(perturbedModels{i}).comparisonStats.accuracy(2)), ...
           round(100 * perturbedModelsComparison.(perturbedModels{i}).comparisonStats.accuracy(3)), ...
           round(perturbedModelsComparison.(perturbedModels{i}).comparisonStats.wEuclidNorm(1)), ...
           round(perturbedModelsComparison.(perturbedModels{i}).comparisonStats.wEuclidNorm(2)), ...
           round(perturbedModelsComparison.(perturbedModels{i}).comparisonStats.wEuclidNorm(3))];
    end
    
    % Table
    iDopaNeuro1CoreValidationTable = table(columnsTable(:, 1), ...
        columnsTable(:, 2), columnsTable(:, 3), ...
        ...
        'VariableNames', ...
        perturbedModels(1:3),...
        'RowNames',...
        {'qualitativeBoth (%)'; ...
        'qualitativeModelSec (%)';...
        'qualitativeModelUpt (%)';...
        'quantitativeBoth (μmol/gDW/hr)';...
        'quantitativeModelSec (μmol/gDW/hr)';...
        'quantitativeModelUpt (μmol/gDW/hr)'});
    display(iDopaNeuro1CoreValidationTable)
    
    
end


