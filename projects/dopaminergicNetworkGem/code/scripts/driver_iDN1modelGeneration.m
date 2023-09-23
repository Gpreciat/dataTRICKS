% Generate iDopaNeuro model varations with different approaches
% Note: It is important to check if the outputDir is correct
% created models e.g.: ~/work/sbgCloud/programExperimental/projects/xomics/results/iDN1

clear

[solverOK, solverInstalled] = changeCobraSolver('gurobi','all');

% Define results directory
if ~contains(char(java.lang.System.getProperty('user.name')),'rfleming')
outputDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1_final'];
else
   outputDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'iDopaNeuro1']; 
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

% Load Recon3D (noLumpedRxns)
inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
    filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input.mat';
load([inputFolder filesep genericModelName])

% SPECIFIC DATA
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
% bibliomicData
bibliomicData = 'bibliomicData.xlsx';
specificData = preprocessingOmicsModel([dataFolder bibliomicData], 1, 1);
% exoMet
specificData.exoMet = readtable([dataFolder 'exoMet']);
% transcriptomicData
specificData.transcriptomicData = readtable([dataFolder 'transcriptomicData']);
specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);

% Parameters tested
% "thermoKernel_oneRxnPerActiveGene_transcriptomicsT2_limitBoundary.100000_inactiveGenesT_openIons_omicsOverCuration"
param.tissueSpecificSolver = 'fastCore';
param.activeGenesApproach = 'oneRxnPerActiveGene'; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
param.transcriptomicThreshold = 0; % [0 1 2];
param.TolMinBoundary = -1e5;
param.TolMaxBoundary = 1e5;
param.inactiveGenesTranscriptomics = false; % [true false];
param.closeIons = true; % [true false];
param.curationOverOmics = false;

% General parameters
param.workingDirectory = outputDir;
param.printLevel = 2;
param.setObjective = ''; % No objective function
feasTol = getCobraSolverParams('LP', 'feasTol');
param.boundPrecisionLimit = feasTol * 10;
param.fluxEpsilon = feasTol * 10;
param.fluxCCmethod = 'fastcc';
param.weightsFromOmics = 1;
param.metabolomicWeights = 'mean';
%param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true; % Cell culture information
param.debug = true;

% Start the diary
if isunix()
    name = getenv('USER');
else
    name = getenv('username');
end
param.diaryFilename = [outputDir filesep datestr(now,30) '_' name '_diary.txt'];

% Remove expressionRxns
if isfield(model, 'expressionRxns')
    model = rmfield(model, 'expressionRxns');
end

% Create models
[iDopaNeuro1, ~] = XomicsToModel(model, specificData, param);

save([outputDir filesep 'iDopaNeuro1.mat'], 'iDopaNeuro1')
