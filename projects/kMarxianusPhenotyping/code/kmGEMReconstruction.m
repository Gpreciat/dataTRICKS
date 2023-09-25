% MATLAB Script: Kluyveromyces marxianus reconstruction

% Description:
% This script aims to improve the genome-scale metabolic model of the bacterium
% Kluyveromyces marxianus by addressing connectivity gaps and achieving
% stoichiometric, flux, and thermodynamic consistency. The initial model was
% generated from the genetic sequence of Kluyveromyces marxianus, but it lacks
% proper connectivity. To enhance the model, we will integrate information from
% two previously published models of Kluyveromyces marxianus and Saccharomyces cerevisiae.
% By combining these models, we intend to reconstruct a consistent metabolic
% network, filling gaps, and resolving inconsistencies. The improved model
% will be used for predictions under various stress conditions to verify its
% accuracy and biological relevance.

% Steps:
% 1. Load Kluyveromyces marxianus models.
% 2. Identifying Common Reactions and Metabolites.


% 2. Generate a genome-scale metabolic model for Kluyveromyces marxianus (Model KM).
% 3. Identify and quantify connectivity gaps in Model KM.
% 4. Retrieve relevant metabolic information from the previously published models
%    of Kluyveromyces marxianus and Saccharomyces cerevisiae.
% 5. Combine the information from the three models to create an integrated metabolic
%    model for Kluyveromyces marxianus with improved connectivity.
% 6. Validate the improved model for stoichiometric, flux, and thermodynamic consistency.
% 7. Simulate and predict metabolic behaviors under different stress conditions.
% 8. Perform comparisons with experimental data to assess model accuracy and performance.

% Author: German Preciat
% Date: 05/30/2023

% Note: This script requires the Cobratoolbox and other relevant functions to be installed.
%       Make sure to have the necessary dependencies and data files available.

%% Step 1: Load Kluyveromyces marxianus models.
% Loads the metabolic models of Kluyveromyces marxianus from two different file formats:
% 1. XML file containing the model with gaps, kmUdg.
% 2. MAT file containing the base model, kmGEMv1.

clear

% Define directories
filePath = regexprep(matlab.desktop.editor.getActiveFilename, ['code' filesep 'kmGEMReconstruction.m'], '');
dataDir = [filePath 'data' filesep];  
resultsDir = [filePath 'results' filesep];  

if ~isfile([resultsDir '1modelsLoaded.mat'])
    
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
    kmUdg.rxnFormulas = printRxnFormula(kmUdg, 'printFlag', 0);
    
    % kmGEMv1 model
    % Size of S in kmGEMv1 model: 1531 x 1913
    load([dataDir 'km.mat'])
    kmGEMv1 = model;
    kmGEMv1.rxnFormulas = printRxnFormula(kmGEMv1, 'printFlag', 0);
    
    % Save debug file
    clear model
    save([resultsDir '1modelsLoaded.mat'])
    
end

%% Step 2: Identifying Common Reactions and Metabolites
% Identifies the common reactions and metabolites between the two Kluyveromyces 
% marxianus models using available identifiers for further analysis.

clearvars -except resultsDir

if ~isfile([resultsDir '2ndSection.mat'])
    
    load([resultsDir '1modelsLoaded.mat'])
    
    [rxnTable, metTable]  = commonIdsBetweenModels(kmUdg, kmGEMv1, 'both', true);
    
    nRows = length(metsInCommon) * length(objectives);
varTypes = {'double', 'string', 'string', 'string', 'string', 'string',...
    'string', 'string', 'string', 'string', 'double', 'double', 'double',...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'logical', 'string'};
varNames = {'modelId', 'dirName', 'tissueSpecificSolver', 'activeGenesApproach',...
    'transcriptomicThreshold', 'limitBounds', 'genesTranscriptomics', 'ions', ...
    'preferenceCurationOrOmics', 'objective', 'qualitativeBoth', 'quantitativeBoth',...
    'spearmanBoth', 'qualitativeModelSec', 'quantitativeModelSec', 'spearmanModelSec', ...
    'quantitativeModelUpt', 'qualitativeModelUpt', 'spearmanModelUpt', 'nOfmets', ...
    'nOfrxns', 'rankOfS', 'ATPtm', 'ME'};
multidimensionalComparisonStats = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
    'VariableNames', varNames);
    
    % Add fields and IDs from the original model
    sameKeggRxnId = kmUdg.rxns(ismember(kmUdg.rxnKEGGID, kmGEMv1.rxnKEGGID));
    sameKeggRxnId2 = kmGEMv1.rxns(ismember(kmGEMv1.rxnKEGGID, kmUdg.rxnKEGGID));
    [kmGEMv1.rxnBiGGID, kmGEMv1.rxnBioCycID, kmGEMv1.rxnMetaNetXID, kmGEMv1.rxnRheaID] = deal(cell(size(kmGEMv1.rxns)));
    % Save ids in km2
    for i = 1:length(kmUdg.rxnKEGGID)
        ids = split(kmUdg.rxnKEGGID{i}, '; ');
        if length(ids) > 1
            for j = 1:length(ids)
                idBool = ismember(kmGEMv1.rxnKEGGID, ids{j});
                if any(idBool)
                    kmGEMv1.rxnBiGGID{idBool} = kmUdg.rxnBiGGID{i};
                    kmGEMv1.rxnBioCycID{idBool} = kmUdg.rxnBioCycID{i};
                    kmGEMv1.rxnMetaNetXID{idBool} = kmUdg.rxnMetaNetXID{i};
                    kmGEMv1.rxnRheaID{idBool} = kmUdg.rxnRheaID{i};
                end
            end
        else
            idBool = ismember(kmGEMv1.rxnKEGGID, ids);
            if any(idBool) && ~isempty(char(ids))
                kmGEMv1.rxnBiGGID(idBool) = kmUdg.rxnBiGGID(i);
                kmGEMv1.rxnBioCycID(idBool) = kmUdg.rxnBioCycID(i);
                kmGEMv1.rxnMetaNetXID(idBool) = kmUdg.rxnMetaNetXID(i);
                kmGEMv1.rxnRheaID(idBool) = kmUdg.rxnRheaID(i);
            end
        end
    end
    
    % Add fields and IDs from the original model
    [kmGEMv1.metBiGGID, kmGEMv1.metBioCycID, kmGEMv1.metHMDBID, kmGEMv1.metMetaNetXID, ...
        kmGEMv1.metSEEDID, kmGEMv1.metisinchikeyID] = deal(cell(size(kmGEMv1.mets)));
    kmGEMv1.metCharges = zeros(size(kmGEMv1.mets));
    % Save ids in km2
    for i = 1:length(kmUdg.metKEGGID)
        ids = split(kmUdg.metKEGGID{i}, '; ');
        if length(ids) > 1
            for j = 1:length(ids)
                idBool = ismember(kmGEMv1.metKEGGID, ids{j});
                if any(idBool)                    
                    kmGEMv1.metBiGGID(idBool) = kmUdg.metBiGGID(i);
                    kmGEMv1.metBioCycID(idBool) = kmUdg.metBioCycID(i);
                    kmGEMv1.metCharges(idBool) = kmUdg.metCharges(i);
                    kmGEMv1.metHMDBID(idBool) = kmUdg.metHMDBID(i);
                    kmGEMv1.metMetaNetXID(idBool) = kmUdg.metMetaNetXID(i);
                    kmGEMv1.metSEEDID(idBool) = kmUdg.metSEEDID(i);
                    kmGEMv1.metisinchikeyID(idBool) = kmUdg.metisinchikeyID(i);
                end
            end
        else
            idBool = ismember(kmGEMv1.rxnKEGGID, ids);
            if any(idBool) && ~isempty(char(ids))
                kmGEMv1.metBiGGID(idBool) = kmUdg.metBiGGID(i);
                kmGEMv1.metBioCycID(idBool) = kmUdg.metBioCycID(i);
                kmGEMv1.metCharges(idBool) = kmUdg.metCharges(i);
                kmGEMv1.metHMDBID(idBool) = kmUdg.metHMDBID(i);
                kmGEMv1.metMetaNetXID(idBool) = kmUdg.metMetaNetXID(i);
                kmGEMv1.metSEEDID(idBool) = kmUdg.metSEEDID(i);
                kmGEMv1.metisinchikeyID(idBool) = kmUdg.metisinchikeyID(i);
            end
        end
    end
    
    [kmUdg, ~, ~] = removeRxns(kmUdg, sameKeggRxnId);
    % km_ha = alphabetizeModel(km_ha);
    
    %     for i = 1:length(km2.mets)
    %         if ~isempty(km2.mets{i})
    %             km2.mets{i} = km2.metKEGGID{i};
    %         end
    %     end
    %     for i = 1:length(km_ha.mets)
    %         if ~isempty(km_ha.mets{i})
    %             km_ha.mets{i} = km_ha.metKEGGID{i};
    %         end
    %     end
    kmUdg.rxnsFormulas = printRxnFormula(kmUdg, 'printFlag', 0);
    kmGEMv1.rxnsFormulas = printRxnFormula(kmGEMv1, 'printFlag', 0);
    
    % By rxnFormula
    % Change met Id to Keeg ID
    temp_km_ha = kmUdg;
    for i = 1:length(temp_km_ha.metKEGGID)
        if ~isempty(temp_km_ha.metKEGGID{i})
            ids = split(temp_km_ha.metKEGGID{i}, '; ');
            if length(ids) > 1
                for j = 1:length(ids)
                    if ismember(ids{j}, kmGEMv1.metKEGGID)
                        temp_km_ha.mets{i} = ids{j};
                        continue
                    end
                end
            end
        end
    end
    temp_km2 = kmGEMv1;
    temp_km2.mets = temp_km2.metKEGGID;
    temp_km_ha.rxnsFormulas = printRxnFormula(temp_km_ha, 'printFlag', 0);
    temp_km2.rxnsFormulas = printRxnFormula(temp_km2, 'printFlag', 0);
    
    % Compare formulas
    rxns2Remove = [];
    inKmFormulas = [];
    for i = 1:length(kmUdg.rxns)
        if ismember(temp_km_ha.rxnsFormulas{i}, temp_km2.rxnsFormulas)
            rxns2Remove = [rxns2Remove; kmUdg.rxns(i)];
            inKmFormulas = [inKmFormulas kmGEMv1.rxns(ismember(temp_km2.rxnsFormulas, temp_km_ha.rxnsFormulas{i}))'];
        end
    end
    if ~isempty(rxns2Remove)
        [kmUdg, ~, ~] = removeRxns(kmUdg, rxns2Remove);
    end
    
    % By EC number (none)
    rxns2Remove = [];
    EcNumberRepeated = [];
    ecnumber = [];
    for i = 1:length(kmGEMv1.rxnECNumbers)
        rxnECNumberskm2 = strtrim(split(kmGEMv1.rxnECNumbers{i}, ';'));
        for j = 1:length(kmUdg.rxnECNumbers)
            rxnECNumbersKm_ha = strtrim(split(kmUdg.rxnECNumbers{j}, '; '));
            for k = 1:length(rxnECNumberskm2)
                if ismember(rxnECNumberskm2(k), rxnECNumbersKm_ha)  && ~isempty(rxnECNumberskm2{k})
                    rxns2Remove = [rxns2Remove kmUdg.rxns(j)];
                    EcNumberRepeated = [EcNumberRepeated; kmGEMv1.rxns(i)];
                end
            end
        end
    end
    if ~isempty(rxns2Remove)
        [kmUdg, ~, ~] = removeRxns(kmUdg, rxns2Remove);
    end
    
    % Add prefix to reactions in both
    reactionsInBoth = findRxnIDs(kmGEMv1, unique([sameKeggRxnId2; EcNumberRepeated; inKmFormulas]));
    for i = 1:length(reactionsInBoth)
        kmGEMv1.rxns{reactionsInBoth(i)} = ['both_' kmGEMv1.rxns{reactionsInBoth(i)}];
    end
    
    % Save debug file
    save([resultsDir '2ndSection.mat'])
    
end

%%

clearvars -except outputDir

if ~isfile([resultsDir '3rdSection.mat'])
    
    clearvars -except outputDir
    
    load([resultsDir '2ndSection.mat'])
        
    kmUdg.mets = regexprep(kmUdg.mets, '\[c0\]', '\[c\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[e0\]', '\[e\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[z5\]', '\[n\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[z6\]', '\[p\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[g0\]', '\[g\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[m0\]', '\[m\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[r0\]', '\[r\]');
    kmUdg.mets = regexprep(kmUdg.mets, '\[z7\]', '\[v\]');
    kmUdg.rxnsFormulas = printRxnFormula(kmUdg, 'printFlag', 0);
    
    for i = 1:length(kmUdg.rxns)
        if ~ismember(['origKmha_' kmUdg.rxns{i}], kmGEMv1.rxns)
            kmGEMv1 = addReaction(kmGEMv1, ['origKmha_' kmUdg.rxns{i}], ...
                'reactionFormula', kmUdg.rxnsFormulas{i}, ...
                'reactionName', kmUdg.rxnNames{i}, ...
                'lowerBound', kmUdg.lb(i), ...
                'upperBound', kmUdg.ub(i), ...
                'geneRule', kmUdg.grRules{i});
            idx = length(kmGEMv1.rxns);
            kmGEMv1.rxnECNumbers(idx) = kmUdg.rxnECNumbers(i);
            kmGEMv1.rxnBiGGID(idx) = kmUdg.rxnBiGGID(i);
            kmGEMv1.rxnBioCycID(idx) = kmUdg.rxnBioCycID(i);
            kmGEMv1.rxnMetaNetXID(idx) = kmUdg.rxnMetaNetXID(i);
            kmGEMv1.rxnRheaID(idx) = kmUdg.rxnRheaID(i);
        end
    end
    
    kmGEMv1 = alphabetizeModel(kmGEMv1);
    kmGEMv1.rxnsFormulas = printRxnFormula(kmGEMv1, 'printFlag', 0);
    
    % Save debug file
    save([resultsDir '3rdSection.mat'])
    
end
%%

load([resultsDir '3rdSection.mat'])


    
    % Cambiar solver
    [~, ~] = changeCobraSolver('mosek', 'all');
specificData.activeReactions = kmGEMv1.rxns(contains(kmGEMv1.rxns, 'origKmha_') | contains(kmGEMv1.rxns, 'both_'));
% specificData.transcriptomicData = table;
% specificData.transcriptomicData.genes = km2.genes(randi(length(km2.genes), [30, 1]));
% specificData.transcriptomicData.expVal = 200 * rand(30, 1);

% XomicsToModel
% General parameters
param.workingDirectory = resultsDir;
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
[km_haNew, ~] = XomicsToModel(kmGEMv1, specificData, param);

% Salvar modelo
save([resultsDir filesep 'km_haNew.mat'], 'km_haNew')

%%

% verificar xilosa

km_haNew = changeObjective(km_haNew, km_haNew.rxns{1});
fbaKm = optimizeCbModel(km_haNew, 'max');
