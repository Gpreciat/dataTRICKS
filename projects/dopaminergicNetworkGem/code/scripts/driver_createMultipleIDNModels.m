% Generate iDopaNeuro model conditions with different conditions using the
% function XomicsToMultipleModels. The current conditions are divided in 9
% groups:
%
%   1. Generic model: Recon3DModel_301_xomics_input
%   2. specificData: ATPMlb186 ATPMlb372 ATPMlb744 ATPMlb186ATPtm ATPMlb372ATPtm ATPMlb744ATPtm
%   3. cobraSolver: gurobi
%   4. tissueSpecificSolver: fastCore - thermoKernel
%   5. activeGenesApproach: deleteModelGenes - oneRxnPerActiveGene
%   6. transcriptomicThreshold: [0 2]
%   7. limitBounds: [1e3 1e4]
%   8. inactiveGenesTranscriptomics: [true false]
%   9. closeIons: [true false]
%  10. curationOverOmics: [true false]
%
% The models are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1
%% Set directory and specificData

clear

% specificData
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
bibliomicData = 'bibliomicData.xlsx';
specificData = preprocessingOmicsModel([dataFolder bibliomicData], 1, 1);
specificData.exoMet = readtable([dataFolder 'exoMet']);
specificData.transcriptomicData = readtable([dataFolder 'transcriptomicData']);
specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);

%% Multiple options

% Generic model
inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
    filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input.mat';
load([inputFolder filesep genericModelName])
modelGenerationConditions.genericModel.model = model;

% Data to compare
if ~contains(char(java.lang.System.getProperty('user.name')),'rfleming')

    modelsDir = '/Users/gpreciat/Desktop/tmp/multidimensionalModelGeneration';
    if ~isfolder(modelsDir)
        mkdir(modelsDir);
    end
    cd(modelsDir)
    
    % Define results directory
    modelGenerationConditions.outputDir = modelsDir;
    % modelGenerationConditions
    % Note: if the data does not vary it is enough to declare them in param.
    % Here the non-varying conditions are left for demonstrative purposes
    modelGenerationConditions.cobraSolver = {'gurobi'}; %  {'ibm_cplex', 'gurobi'};
    modelGenerationConditions.activeGenesApproach = {'deleteModelGenes', 'oneRxnPerActiveGene'}; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
    modelGenerationConditions.transcriptomicThreshold = [0 2]; % [0 1 2];
    modelGenerationConditions.boundsToRelaxExoMet = {'b'}; % {'l', 'b', 'u'}
    modelGenerationConditions.closeIons = [true false]; % [true false];
    modelGenerationConditions.tissueSpecificSolver = {'fastCore', 'thermoKernel'}; % {'fastCore', 'thermoKernel'};
    modelGenerationConditions.inactiveGenesTranscriptomics = [true false]; % [true false];
    modelGenerationConditions.limitBounds = [1e3 1e4 1e5]; % [1e3 1e4 1e5];
    modelGenerationConditions.curationOverOmics = [true false]; % [true false];
    
    % Different lb for the demand of ATP (ATPM)
    ATPMidx = ismember(specificData.rxns2constrain.rxnID, 'ATPM');
    modelGenerationConditions.specificData.ATPMlb186 = specificData;
    modelGenerationConditions.specificData.ATPMlb372 = specificData;
    modelGenerationConditions.specificData.ATPMlb372.rxns2constrain.lb(ATPMidx) = 372;
    modelGenerationConditions.specificData.ATPMlb744 = specificData;
    modelGenerationConditions.specificData.ATPMlb744.rxns2constrain.lb(ATPMidx) = 744;
    
    % Force activity of ATPtm (transport of ATP from [m] to [c])
    ATPtmidx = length(specificData.rxns2constrain.rxnID) + 1;
    specificData2 = specificData;
    specificData2.rxns2constrain.rxnID(ATPtmidx) = {'ATPtm'};
    specificData2.rxns2constrain.lb(ATPtmidx) = 0.0001;
    specificData2.rxns2constrain.ub(ATPtmidx) = 10000;
    specificData2.rxns2constrain.constraintDescription(ATPtmidx) = {'Force activity'};
    modelGenerationConditions.specificData.ATPMlb186ATPtm = specificData2;
    modelGenerationConditions.specificData.ATPMlb372ATPtm = specificData2;
    modelGenerationConditions.specificData.ATPMlb372ATPtm.rxns2constrain.lb(ATPMidx) = 372;
    modelGenerationConditions.specificData.ATPMlb744ATPtm = specificData2;
    modelGenerationConditions.specificData.ATPMlb744ATPtm.rxns2constrain.lb(ATPMidx) = 744;
    
else
    modelsDir = '~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/iDopaNeuro1b';
    if ~isfolder(modelsDir)
        mkdir(modelsDir);
    end
    cd(modelsDir)

    % Define results directory
    modelGenerationConditions.outputDir = modelsDir;
    
    %debug generation of a single model type
    % modelGenerationConditions
    % Note: if the data does not vary it is enough to declare them in param.
    % Here the non-varying conditions are left for demonstrative purposes
    modelGenerationConditions.cobraSolver = {'gurobi'}; %  {'ibm_cplex', 'gurobi'};
    modelGenerationConditions.activeGenesApproach = {'oneRxnPerActiveGene'}; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
    modelGenerationConditions.transcriptomicThreshold = 0; % [0 1 2];
    modelGenerationConditions.boundsToRelaxExoMet = {'b'}; % {'l', 'b', 'u'}
    modelGenerationConditions.closeIons = false; % [true false];
    modelGenerationConditions.tissueSpecificSolver = {'thermoKernel'}; % {'fastCore', 'thermoKernel'};
    modelGenerationConditions.inactiveGenesTranscriptomics = true; % [true false];
    modelGenerationConditions.limitBounds = 1000; % [1e3 1e4 1e5]; %was inf
    modelGenerationConditions.curationOverOmics = false; % [true false];
        
    param.printLevel = 2;
end

% specificData
modelGenerationConditions.specificData.specificData = specificData;

%% Fixed options

param.setObjective = ''; % No objective function
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';
param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true; % Cell culture information
param.debug = false;

%% Create models

% Remove expressionRxns
if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end

% Create models
replaceModels = false;
directoriesWithModels = XomicsToMultipleModels(modelGenerationConditions, param);

display('Models saved in:')
display(modelsDir)