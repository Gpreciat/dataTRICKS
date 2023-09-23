% Driver to test the objective functions unWeightedTCBMflux and 
% unWeightedTCBMfluxConc using the models saved in:
%
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/multidimensionalModelGeneration

clear

% Define results directory
modelsDir = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'multidimensionalModelGeneration'];

% Get a list of all files and folders in this folder.
directoriesWithModels = dir(modelsDir);
directoriesWithModels = struct2cell(directoriesWithModels([directoriesWithModels.isdir]));
directoriesWithModels = directoriesWithModels(1, 3:end)';

%% Flux estimation (unWeightedTCBMflux and unWeightedTCBMfluxConc)

% Select solver
[~, ~] = changeCobraSolver('gurobi', 'LP');

c = 0;
for i = 1:size(directoriesWithModels, 1)
    i = 1;
    workingDirectory = [modelsDir filesep directoriesWithModels{i} filesep];
    
    % Identify number dimensions
    load([workingDirectory 'Model.mat'])
    
	% unWeightedTCBMflux
    model = Model;
    model.osenseStr = 'min';
    model.cf = 0;
    model.cr = 0;
    model.g = 2;
    model.u0 = 0;
    model.f = 1;
    tcbmParam.method = 'fluxes';
    tcbmParam.printLevel = 1;
    tcbmParam.solver = 'mosek';
    [unWeightedTCBMflux, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);

    % unWeightedTCBMfluxConc
    model = Model;
    model.osenseStr = 'min';
    model.cf = 0;
    model.cr = 0;
    model.g = 2;
    model.u0 = 0;
    model.f = 1;
    tcbmParam.method = 'fluxConc';
    tcbmParam.printLevel = 1;
    tcbmParam.solver = 'mosek';
    [unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);

end
