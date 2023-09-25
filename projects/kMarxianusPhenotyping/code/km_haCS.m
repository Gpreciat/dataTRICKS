% Driver to create a context-specific genome-scale model of km_ha cell line
% based on transcriptomic data
%
% Meast model
% https://github.com/SysBioChalmers/yeast-GEM
%
% Sources of interest:
% https://www.nature.com/articles/s41467-021-25158-6#data-availability

%% Initialise

clear

% Define directories
filePath = regexprep(matlab.desktop.editor.getActiveFilename, 'km_haCS.m', '');
dataDir = [filePath filesep 'data' filesep];
outputDir = [filePath 'results' filesep];

if ~isfile([outputDir '1stSection.mat'])
    
    % Read models
    
    % km_ha
    km_ha = readCbModel([dataDir 'km_ha.xml'],'fileType','SBML','defaultBound', 1000);
    km_ha.grRules = km_ha.rules;
    % yeast model
    load([dataDir 'yeast-GEM.mat'])
    yeastModel = model;
    % Size of S in km_ha:       1228 x 1053
    % Size of S in yeast model: 2744 x 4063
    
    % Cambiar solver
    [~, ~] = changeCobraSolver('mosek', 'all');
    
    clear model
    save([outputDir '1stSection.mat'])
    
end
%% Consistent naming of mets and rxns in both models (BiGG ids)

clearvars -except outputDir
load([outputDir '1stSection.mat'])

% Add BiGG ids in the yeast model
% Read Ids
yeastBiggRxns = readtable([dataDir 'BiGGrxnDictionary.csv']);
yeastBiggMets = readtable([dataDir 'BiGGmetDictionary.csv']);
% mets
[~, locb] = ismember(yeastModel.mets, regexprep(yeastBiggMets.mets, '(\[\w\])', ''));
yeastModel.metBiGGID = cell(length(locb), 1);
for i = 1:length(locb)
    yeastModel.metBiGGID{i} = '';
    if locb(i) ~= 0
        yeastModel.metBiGGID(i) = yeastBiggMets.bigg(locb(i));
    end
end
% rxns
[~, locb] = ismember(yeastModel.rxns, yeastBiggRxns.rxns);
yeastModel.rxnBiGGID = cell(length(locb), 1);
for i = 1:length(locb)
    yeastModel.rxnBiGGID{i} = '';
    if locb(i) ~= 0
        yeastModel.rxnBiGGID(i) = yeastBiggRxns.bigg(locb(i));
    end
end
% Change compartments in the yeast model
% km_ha compartmets:
%    {'e0'}
%    {'c0'}
%    {'m0'}
%    {'z6'}
%    {'z5'}
%    {'r0'}
%    {'z7'}
%    {'g0'}
yeastModel.comps(ismember(yeastModel.comps, [{'c'}; {'ce'}])) = {'c0'};
yeastModel.comps(ismember(yeastModel.comps, {'e'})) = {'e0'};
yeastModel.comps(ismember(yeastModel.comps, [{'m'}; {'mm'}])) = {'m0'};
yeastModel.comps(ismember(yeastModel.comps, [{'g'}; {'gm'}])) = {'g0'};
yeastModel.comps(ismember(yeastModel.comps, [{'er'}; {'erm'}])) = {'r0'};

% Change mets ids to BiGG id and add the correct compartment in both models
% Yeast
for i = 1:length(yeastModel.mets)
    if ~isempty(yeastModel.metBiGGID{i})
        yeastModel.mets{i} = [yeastModel.metBiGGID{i} '[' yeastModel.comps{yeastModel.metComps(i)} ']'];
    else
        yeastModel.mets{i} = [yeastModel.mets{i} '[' yeastModel.comps{yeastModel.metComps(i)} ']'];
    end
end
% km_ha
for i = 1:length(km_ha.mets)
    comp = regexp(km_ha.mets{i}, '(\[\w\d\])', 'match', 'once');
    if ~isempty(km_ha.metBiGGID{i})
        mets = split(km_ha.metBiGGID{i}, '; ');
        met = mets{1};
        if length(mets) > 1
            metIdx = find(ismember(mets, yeastModel.metBiGGID));
            if ~isempty(metIdx)
                met = mets{metIdx};
            end
        end
        km_ha.mets{i} = [met comp];
    end
end

% Change rxns ids to BiGG id and add the correct compartment in both models
% Yeast
for i = 1:length(yeastModel.rxns)
    if ~isempty(yeastModel.rxnBiGGID{i})
        yeastModel.rxns{i} = yeastModel.rxnBiGGID{i};
    end
end
% Create a new model with joint rxns
model = createModel();
for i = 1:length(km_ha.rxns)
    if ~isempty(km_ha.rxnBiGGID{i})
        rxns = split(km_ha.rxnBiGGID{i}, '; ');
        for j = 1:length(rxns)
            idx = findRxnIDs(yeastModel, rxns{j})
        
        
        rxn = rxns{1};
        if length(rxns) > 1
            rxnIdx = find(ismember(rxns, yeastModel.rxnBiGGID));
            if ~isempty(rxnIdx)
                for j = 1:length(rxnIdx)
                
                    model = addReaction(model, yeastModel.rxns{rxnIdx(i)}, 'reactionFormula', yeastModel.rxnsFormulas{i});
                    km_ha.lb(end) = yeastModel.lb(i);
                    km_ha.ub(end) = yeastModel.ub(i);
                    km_ha.rxnECNumbers(end) = yeastModel.eccodes(i);
                    km_ha.grRules(end) = yeastModel.grRules(i);
                    km_ha.rxnBiGGID(end) = yeastModel.rxnBiGGID(i);
                    
                    
                    
                end
                rxn = rxns{rxnIdx};
            end
        end
        km_ha.rxns{i} = rxn;
    end
end

%%

yeastModel.rxnsFormulas = printRxnFormula(yeastModel, 'printFlag', 0);
km_ha.rxnsFormulas = printRxnFormula(km_ha, 'printFlag', 0);

% Rename metabolites of the yeas model based on the EC numbers
for i = 1:length(km_ha.rxnECNumbers)
    % Split joint EC numbers and select one
    ecNumbers = split(km_ha.rxnECNumbers{i}, '; ');
    if length(ecNumbers) > 1
        for j = 1:length(ecNumbers)
            if ismember(ecNumbers(j), yeastModel.eccodes)
                km_ha.rxnECNumbers{i} = ecNumbers(j);
            end
        end
    end
    % Compare rxnFormulas for same EC number
    [bool locb] = ismember(yeastModel.eccodes, km_ha.rxnECNumbers{i});
    if sum(bool) > 1 && ~isempty(km_ha.rxnECNumbers{i})
        
        km_ha.rxnsFormulas{i}
        yeastModel.rxnsFormulas(find(bool))
        display('aqui')
    end
end
    
   if ecSameBool(i)
       isequal(km_ha.rxnECNumbers, yeastModel.eccodes(i))
   end


regexprep(yeastBiggMets.mets, '(\[\w\])', '')

%% Mix models

%[SConsistentMetBool1, SConsistentRxnBool1, SInConsistentMetBool1, ...
%   SInConsistentRxnBool1, unknownSConsistencyMetBool1, ...
%    unknownSConsistencyRxnBool1, model1, stoichConsistModel1] = ...
%    findStoichConsistentSubset(km_ha);
%
%display(['SConsistentMetBool1 ' num2str(sum(SConsistentMetBool1))])

bool = ~ismember(yeastModel.eccodes, km_ha.rxnECNumbers);
for i = 1:length(bool)
    if bool(i)
        km_ha = addReaction(km_ha, ['new_' yeastModel.rxns{i}], 'reactionFormula', yeastModel.rxnsFormulas{i});
        km_ha.lb(end) = yeastModel.lb(i);
        km_ha.ub(end) = yeastModel.ub(i);
        km_ha.rxnECNumbers(end) = yeastModel.eccodes(i);
        km_ha.grRules(end) = yeastModel.grRules(i);
        km_ha.rxnBiGGID(end) = yeastModel.rxnBiGGID(i);
    end
end


%[SConsistentMetBool2, SConsistentRxnBool2, SInConsistentMetBool2, ...
%    SInConsistentRxnBool2, unknownSConsistencyMetBool2, ...
%    unknownSConsistencyRxnBool2, model2, stoichConsistModel2] = ...
%    findStoichConsistentSubset(km_ha);
%
%display(['SConsistentMetBool2 ' num2str(sum(SConsistentMetBool2))])

%%

% Leer datos transcriptomicos (inventados)
specificData.transcriptomicData = readtable([dataDir 'transcriptomicData.txt']);
specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);
specificData.transcriptomicData.genes(1:length(km_ha.genes)) = km_ha.genes;

% XomicsToModel
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
param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
param.addCoupledRxns = 1;
param.nonCoreSinksDemands = 'closeAll';
param.closeUptakes = true; % Cell culture information
param.debug = true;
param.tissueSpecificSolver = 'thermoKernel';
param.activeGenesApproach = 'oneRxnPerActiveGene'; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
param.transcriptomicThreshold = 2; % [0 1 2];
param.TolMinBoundary = -1e4;
param.TolMaxBoundary = 1e4;
param.inactiveGenesTranscriptomics = true; % [true false];
param.closeIons = false; % [true false];

% Correr la funcion
[km_haNew, ~] = XomicsToModel(km_ha, specificData, param);

% Salvar modelo
save([outputDir filesep 'km_haNew.mat'], 'km_haNew')