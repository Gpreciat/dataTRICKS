% combine_genomic_models.m
%
% Description:
% This MATLAB script aims to combine two genome-scale metabolic models of
% Kluyveromyces marxianus (K. marxianus) obtained from different sources:
% one generated using the MERLIN software with sequencing data and another
% downloaded from a literature publication titled "Reconstruction and
% analysis of a Kluyveromyces marxianus genome-scale metabolic model"
% (published in 2019). The goal is to create a unified and comprehensive
% metabolic model of K. marxianus by integrating the information from both
% models while avoiding redundancy.
%
% The first model (Model 1) from MERLIN contains multiple identifiers for
% metabolites and reactions, as well as several cellular compartments.
% However, it is disconnected and has many gaps in the network. On the
% other hand, the second model (Model 2) from the literature is better
% connected, but it has fewer cellular compartments.
%
% Strategy:
% 1. Load both models, Model 1 (MERLIN-generated) and Model 2
%    (Literature-based).
%
% 2. Merge the cellular compartments from both models, ensuring no
%    duplication of compartments.
%
% 3. Handle redundant and unique metabolite and reaction identifiers:
%    a) Remove duplicate metabolite and reaction identifiers from Model 1
%       that have matches in Model 2.
%    b) Add unique metabolite and reaction identifiers from Model 2 to
%       Model 1.
%
% 4. Address stoichiometry differences for overlapping reactions between
%    the two models. If an overlapping reaction exists with the same EC
%    numbers, update the stoichiometry using information from Model 2.
%
% 5. Combine other relevant information (e.g., gene-protein-reaction
%    associations, reaction directions, etc.) from both models.
%
% 6. Perform semi-automatic gap-filling using data from both models to
%    ensure the model can produce all biomass components under a minimal
%    medium condition.
%
% 7. Perform manual curation for gene associations, EC numbers, metabolite
%    names, and reaction elemental/charge balance to improve the accuracy
%    of the final model.
%
% 8. Validate the merged model using constraint-based flux balance
%    analysis (FBA) and compare the predictions against experimental data
%    to ensure its accuracy.
%
% 9. The final unified model will be saved in a separate file for further
%    analysis and simulations.

% Author: [Your Name]
% Date: [Date of Script Creation]

%% Identificar las reacciones y metabolitos únicos en cada modelo:
% 1. Obtén las listas de reacciones y metabolitos únicos en cada modelo,
% eliminando duplicados y considerando identificadores y nombres
% específicos de cada modelo.
% 2. Crea una lista con las reacciones y metabolitos que no están presentes
% en ambos modelos.

clear

% Define directories
filePath = regexprep(matlab.desktop.editor.getActiveFilename, ['code\' filesep 'mergeKmModels.m'], '');
dataDir = [filePath 'data' filesep];
resultsDir = [filePath 'results' filesep];

% Read models

% kmUdg
% Size of S in kmUdg:       1228 x 1053
kmUdg = readCbModel([dataDir 'km_ha.xml'], 'fileType', 'SBML', 'defaultBound', 1000);
kmUdg.grRules = kmUdg.rules;
% Rename the compartments in the kmUdg model
kmUdg.mets = regexprep(kmUdg.mets, '\[c0\]', '\[c\]'); %
kmUdg.mets = regexprep(kmUdg.mets, '\[e0\]', '\[e\]');
kmUdg.mets = regexprep(kmUdg.mets, '\[z5\]', '\[n\]');
kmUdg.mets = regexprep(kmUdg.mets, '\[z6\]', '\[p\]');
kmUdg.mets = regexprep(kmUdg.mets, '\[g0\]', '\[g\]');
kmUdg.mets = regexprep(kmUdg.mets, '\[m0\]', '\[m\]');
kmUdg.mets = regexprep(kmUdg.mets, '\[r0\]', '\[r\]');
kmUdg.mets = regexprep(kmUdg.mets, '\[z7\]', '\[v\]');
kmUdg = alphabetizeModel(kmUdg);
kmUdg.rxnFormulas = printRxnFormula(kmUdg, 'printFlag', 0);

kmUdgFBA = optimizeCbModel(kmUdg, 'max');

% kmGEMv1 model
% Size of S in kmGEMv1 model: 1531 x 1913
load([dataDir 'km.mat'])
kmGEMv1 = model;
kmGEMv1 = alphabetizeModel(kmGEMv1);
kmGEMv1.rxnFormulas = printRxnFormula(kmGEMv1, 'printFlag', 0);


kmGEMv1FBA = optimizeCbModel(kmGEMv1, 'max');

% yeast model
% % Size of S in kmGEMv1 model: 1531 x 1913
% load([dataDir 'km.mat'])
% kmGEMv1 = model;
% kmGEMv1 = alphabetizeModel(kmGEMv1);
% kmGEMv1.rxnFormulas = printRxnFormula(kmGEMv1, 'printFlag', 0);

%% Table per compartment

comps = unique(regexprep(kmUdg.mets, '.*\[(.*?)\]', '$1'));
nRows = length(comps);
varTypes = {'string', 'double', 'double', 'double', 'double'};
varNames = {'comp', 'metsUdg', 'rxnsUdg', 'metsGem', 'rxnsGem'};
compsTable = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
    'VariableNames', varNames);

c = 0;
for i = 1:length(comps)
    compsTable.comp{i} = comps{i};
    compModel1 = extractCompModel(kmUdg, comps{i});
    compsTable.metsUdg(i) = length(compModel1.mets);
    compsTable.rxnsUdg(i) = length(compModel1.rxns);
    compModel2 = extractCompModel(kmGEMv1, comps{i});
    if ~isempty(compModel2)
        c = c + 1;
        sharedComps{c, 1} = comps{i};
        compsTable.metsGem(i) = length(compModel2.mets);
        compsTable.rxnsGem(i) = length(compModel2.rxns);
    end
end

%     comp    metsUdg    rxnsUdg    metsGem    rxnsGem
%     ____    _______    _______    _______    _______
%
%     "c"       553        548        894        924
%     "e"       164        170        183        190
%     "g"        12          4          0          0
%     "m"       195        130        350        255
%     "n"        31         14          0          0
%     "p"        75         63          0          0
%     "r"        23          8         73         30
%     "v"         8          3          0          0

%% Merge models per compartment

models = struct;
c = 0;
metMapping = [];
for j = 1:length(sharedComps)

    compModel1 = extractCompModel(kmUdg, sharedComps{j});
    compModel2 = extractCompModel(kmGEMv1, sharedComps{j});

    for i = 1:length(compModel2.mets)
        if ~isempty(compModel2.metKEGGID{i})
            keggBool = ismember(compModel1.metKEGGID, compModel2.metKEGGID{i});
        else
            keggBool = zeros(size(compModel1.mets));
        end
        if ~isempty(compModel2.metInChIString{i})
            inchiBool = ismember(compModel1.metInChIString, compModel2.metInChIString{i});
        else
            inchiBool = zeros(size(compModel1.mets));
        end
        if ~isempty(compModel2.metChEBIID{i})
            chebiBool = ismember(compModel1.metChEBIID, compModel2.metChEBIID{i});
        else
            chebiBool = zeros(size(compModel1.mets));
        end
        if ~isempty(compModel2.metNames{i})
            %             namesBool = ismember(compModel1.metNames, compModel2.metNames{i});
            namesBool = zeros(size(compModel1.mets));
        else
            namesBool = zeros(size(compModel1.mets));
        end
        metBool = keggBool | inchiBool | chebiBool | namesBool;
        if length(find(metBool)) == 1
            c = c + 1;
            metMapping{c, 1} = compModel1.mets{metBool};
            metMapping{c, 2} = compModel2.mets{i};
            compModel1.mets{metBool} = compModel2.mets{i};
        elseif length(find(metBool)) > 1
            error('check')
        end
    end
    [modelNew] = mergeTwoModels(compModel1, compModel2);
    modelNew = alphabetizeModel(modelNew);
    modelNew.rxnFormulas = printRxnFormula(modelNew, 'printFlag', 0);

    % Find unique values and their indices
    [uniqueValues, ~, index] = unique(modelNew.rxnFormulas);
    % Find indices of repeated values & Find logical index for repeated values
    rxnsRepeated = modelNew.rxns(ismember(index, find(histcounts(index, numel(uniqueValues)) > 1)));
    rxns2remove = rxnsRepeated(contains(rxnsRepeated, ['_' sharedComps{j} '0']));
    % Keep IDs
    for k = 1:length(rxns2remove)
        idx = findRxnIDs(modelNew, rxns2remove{k});
        idIdx = setdiff(find(ismember(modelNew.rxnFormulas, modelNew.rxnFormulas(idx))), idx);
        modelNew.rxnRheaID(idIdx) = modelNew.rxnRheaID(idx);
        modelNew.rxnMetaNetXID(idIdx) = modelNew.rxnMetaNetXID(idx);
        modelNew.rxnBioCycID(idIdx) = modelNew.rxnBioCycID(idx);
        modelNew.rxnBiGGID(idIdx) = modelNew.rxnBiGGID(idx);
    end
    % Remove repeated Idx
    modelNew = removeRxns(modelNew, rxns2remove);

    models.(sharedComps{j}) = modelNew;
end

%% Add transport reactions

% kmGEMv1
matches = regexp(kmGEMv1.rxnFormulas, '\[(.*?)\]', 'tokens');
rxnToExtract = [];
for i = 1:length(matches)
    rxnComp = [matches{i}];
    rxnComp = unique([rxnComp{:}]);
    if length(rxnComp) > 1
        rxnToExtract = [rxnToExtract; kmGEMv1.rxns(i)];
    end
end
transportGEM = extractSubNetwork(kmGEMv1, rxnToExtract);

% kmUdg
matches = regexp(kmUdg.rxnFormulas, '\[(.*?)\]', 'tokens');
rxnToExtract = [];
for i = 1:length(matches)
    rxnComp = [matches{i}];
    rxnComp = unique([rxnComp{:}]);
    if length(rxnComp) > 1
        rxnToExtract = [rxnToExtract; kmUdg.rxns(i)];
    end
end
transportUdg = extractSubNetwork(kmUdg, rxnToExtract);

for i = 1:length(metMapping)
    if findMetIDs(transportUdg, metMapping{i, 1}) ~= 0
        transportUdg.mets(findMetIDs(transportUdg, metMapping{i, 1})) = metMapping(i, 2);
    end
end

% Comps
compsUdG = regexprep(transportUdg.mets, '.*\[(.*?)\]', '$1');

% Map met from different compartment
rareCompsBool = ~ismember(compsUdG, sharedComps);
cCompBool = ismember(regexprep(metMapping(:, 1), '.*\[(.*?)\]', '$1'), {'c'});
for i = 1:length(rareCompsBool)
    if rareCompsBool(i)
        idx = ismember(regexprep(metMapping(:, 1), '(\[\w\])', ''), regexprep(transportUdg.mets{i}, '(\[\w\])', '')) & cCompBool;
        if any(idx)
            transportUdg.mets{i} = [regexprep(metMapping{idx, 2}, '(\[\w\])', '') '[' compsUdG{i} ']'];
        end
    end
end

[modelNew] = mergeTwoModels(transportUdg, transportGEM);
modelNew = alphabetizeModel(modelNew);
modelNew.rxnFormulas = printRxnFormula(modelNew, 'printFlag', 0);

% Find unique values and their indices
[uniqueValues, ~, index] = unique(modelNew.rxnFormulas);
% Find indices of repeated values & Find logical index for repeated values
rxnsRepeated = modelNew.rxns(ismember(index, find(histcounts(index, numel(uniqueValues)) > 1)));
rxn2keep = rxnsRepeated(cellfun(@isempty, regexp(rxnsRepeated, ['(_[' strjoin(unique(compsUdG), '') ']0)'], 'match')));
rxns2remove = rxnsRepeated(~cellfun(@isempty, regexp(rxnsRepeated, ['(_[' strjoin(unique(compsUdG), '') ']0)'], 'match')));

% Keep IDs
for i = 1:length(rxns2remove)
    idx = findRxnIDs(modelNew, rxns2remove{i});
    idIdx = setdiff(find(ismember(modelNew.rxnFormulas, modelNew.rxnFormulas(idx))), idx);
    modelNew.rxnRheaID(idIdx) = modelNew.rxnRheaID(idx);
    modelNew.rxnMetaNetXID(idIdx) = modelNew.rxnMetaNetXID(idx);
    modelNew.rxnBioCycID(idIdx) = modelNew.rxnBioCycID(idx);
    modelNew.rxnBiGGID(idIdx) = modelNew.rxnBiGGID(idx);
end
% Remove repeated Idx
modelNew = removeRxns(modelNew, rxns2remove);
models.t = modelNew;

%% Generate models for lonley compartments

compsGem = regexprep(kmGEMv1.mets, '.*\[(.*?)\]', '$1');
for i = 1:length(comps)
    if ~ismember(comps{i}, compsGem)
        compModel = extractCompModel(kmUdg, comps{i});

        for j = 1:length(compModel.mets)
            if ~isempty(compModel.metKEGGID{j})
                keggBool = ismember(kmGEMv1.metKEGGID, compModel.metKEGGID{j});
            else
                keggBool = zeros(size(kmGEMv1.mets));
            end
            if ~isempty(compModel.metInChIString{j})
                inchiBool = ismember(kmGEMv1.metInChIString, compModel.metInChIString{j});
            else
                inchiBool = zeros(size(kmGEMv1.mets));
            end
            if ~isempty(compModel.metChEBIID{j})
                chebiBool = ismember(kmGEMv1.metChEBIID, compModel.metChEBIID{j});
            else
                chebiBool = zeros(size(kmGEMv1.mets));
            end
            if ~isempty(compModel.metNames{j})
                %             namesBool = ismember(compModel1.metNames, compModel2.metNames{i});
                namesBool = zeros(size(kmGEMv1.mets));
            else
                namesBool = zeros(size(kmGEMv1.mets));
            end
            metBool = (keggBool | inchiBool | chebiBool | namesBool) & ismember(compsGem, 'c');
            if length(find(metBool)) > 0
                compModel.mets{j} = [regexprep(kmGEMv1.mets{metBool}, '(\[\w\])', '') '[' comps{i} ']'];
            end
        end
        models.(comps{i}) = compModel;
    end
end

%% Merge models

compModels = fieldnames(models);

finalModel = createModel;
for i = 1:length(compModels)
    finalModel = mergeTwoModels(finalModel, models.(compModels{i}));
end

finalModel = alphabetizeModel(finalModel);
finalModel.rxnFormulas = printRxnFormula(finalModel, 'printFlag', 0);

finalModel.lb(ismember(finalModel.lb, -999999)) = -1000;
finalModel.ub(ismember(finalModel.ub, 999999)) = 1000;

%%
transcriptomicData = readtable([dataDir filesep 'kmExpData.txt']);
[model.expressionRxns, ~] = mapExpressionToReactions(finalModel, transcriptomicData);

%% Check consitency

[~, ~] = changeCobraSolver('mosek', 'all');

feasTol = getCobraSolverParams('LP', 'feasTol');
printLevel = 1;
if 1
    paramFluxConsistency.epsilon = feasTol * 10;
    paramFluxConsistency.method = 'fastcc';
    paramFluxConsistency.printLevel = 1;
else
    paramFluxConsistency.printLevel=1;
    paramFluxConsistency.method = 'null_fastcc';
end

[fluxConsistentMetBool, fluxConsistentRxnBool, fluxInConsistentMetBool, fluxInConsistentRxnBool, ~, fluxConsistModel] =...
    findFluxConsistentSubset(finalModel, paramFluxConsistency);

%
%
% % Leer datos transcriptomicos (inventados)
% specificData.transcriptomicData = transcriptomicData;
%
% % XomicsToModel
% % General parameters
% param.workingDirectory = resultsDir;
% param.printLevel = 2;
% param.setObjective = ''; % No objective function
% feasTol = getCobraSolverParams('LP', 'feasTol');
% param.boundPrecisionLimit = feasTol * 10;
% param.weightsFromOmics = 1;
% param.metabolomicWeights = 'mean';
% param.sinkDMinactive = 1; % Set non-core sinks and demands to inactive
% param.addCoupledRxns = 1;
% param.nonCoreSinksDemands = 'closeAll';
% param.closeUptakes = true; % Cell culture information
% param.debug = true;
% param.tissueSpecificSolver = 'thermoKernel';
% param.activeGenesApproach = 'oneRxnPerActiveGene'; % {'deleteModelGenes', 'oneRxnPerActiveGene'};
% param.transcriptomicThreshold = 2; % [0 1 2];
% param.TolMinBoundary = -1e4;
% param.TolMaxBoundary = 1e4;
% param.inactiveGenesTranscriptomics = true; % [true false];
% param.closeIons = false; % [true false];
%
% % Correr la funcion
% [km_haNew, ~] = XomicsToModel(finalModel, specificData, param);
%
%
[~, ~, ~, ~, ~, ~, ~, stoichConsistModel] = findStoichConsistentSubset(fluxConsistModel, ...
    0, printLevel, [], feasTol * 10);

finalModel = changeObjective(finalModel, 'r_1912');
FBAaerobic = optimizeCbModel(finalModel, 'max');


%%

optionsSampling.samplerName = 'CHRR';
optionsSampling.nStepsPerPoint = 1;
optionsSampling.nPointsReturned = 1;
optionsSampling.toRound = 1;

[P_model, ~] =  sampleCbModel(stoichConsistModel, [], 'CHRR', optionsSampling);

% Depending the size of the polytopes  the optimal sampling its is
% selected as  and . This time, we can avoid the rounding step by inputting
% the rounded polytope from the previous round of sampling (this may take
% some time).
optionsSampling.toRound = 0;
optionsSampling.nPointsReturned = 8 * size(P_model.N, 2);
optionsSampling.nStepsPerPoint = 1000; %size(P_model.N, 2)^2;
%[modelSampling,samples,volume] = sampleCbModel(model, sampleFile, samplerName, options, modelSampling)
[~, samples] =  sampleCbModel(stoichConsistModel, [], optionsSampling.samplerName, optionsSampling, P_model);
[~, samples] =  sampleCbModel(stoichConsistModel, [], optionsSampling.samplerName, optionsSampling);
covS = cov(samples');
numOfExchanges = 20;

%[A_rxns, adjustedA_norms, A_norms] = getTopKNormsAdaptive(A,k,SIntRxnBool)
[S_rxns,S_adjusted_norms, S_norms] = getTopKNormsAdaptive(covS, numOfExchanges, model.SIntRxnBool);
plot_data = [S_adjusted_norms S_norms - S_adjusted_norms];
barh(logPlot_data,'stacked');
xlabel('$log_{10} (Euclidean \: Norm)$', 'interpreter', 'latex')
title('D. Priority metabolites', 'FontSize', 14)
set(gca,'yticklabels',S_rxns);
yticks(1:length(met_names))
