% The iDopaNeuro model is perturbed by the inhibition of mitochondrial
% complexes I and V. In a cell culture with dopaminergic neurons, the
% corresponding perturbations are made by treating them with rotenone or
% oligomycin to mimic the inhibition of mitochondrial complexes I and V,
% respectively. Finally, the in silico predictions are validated against
% the experimental data.
%
% The results are saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/iDopaNeuro1
clear

plotParam.labelType = 'metabolite'; %'platform'
plotParam.saveFigures = 1;

% Select a COBRA model
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep ...
    'iDN1', filesep 'iDopaNeuro1'];


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

load([pathSave filesep 'iDopaNeuro1']);
model = iDopaNeuro1;

%exometabolomic constraints
exoMet = model.XomicsToModelSpecificData.exoMet;

% Identify secretions and uptakes present in the model
exoMetRxnsIdx = ismember(exoMet.rxnID, model.rxns);

% Secretions
rxnTargetSecList = exoMet.rxnID(exoMet.mean > 0 & exoMetRxnsIdx);
rxnTargetSec_mean = exoMet.mean(exoMet.mean > 0 & exoMetRxnsIdx);
rxnTargetSec_SD = exoMet.SD(exoMet.mean > 0 & exoMetRxnsIdx);

% Uptakes
rxnTargetUptList = exoMet.rxnID(exoMet.mean < 0 & exoMetRxnsIdx);
rxnTargetUpt_mean = exoMet.mean(exoMet.mean < 0 & exoMetRxnsIdx);
rxnTargetUpt_SD = exoMet.SD(exoMet.mean < 0 & exoMetRxnsIdx);

% Generate modelUpt
modelUpt = model;
% Open secretion bounds
modelUpt = changeRxnBounds(modelUpt, rxnTargetSecList, 0, 'l');
modelUpt = changeRxnBounds(modelUpt, rxnTargetSecList, max(model.ub), 'u');

fixBounds = 0;
if fixBounds
    % Fix uptake bounds
    meanBound = (model.lb(ismember(model.rxns, rxnTargetUptList)) + model.ub(ismember(model.rxns, rxnTargetUptList))) / 2;
    modelUpt.lb(ismember(model.rxns, rxnTargetUptList)) = meanBound - (model.ub(ismember(model.rxns, rxnTargetUptList)) + model.lb(ismember(model.rxns, rxnTargetUptList))) / 20;
    modelUpt.ub(ismember(model.rxns, rxnTargetUptList)) = meanBound + (model.ub(ismember(model.rxns, rxnTargetUptList)) + model.lb(ismember(model.rxns, rxnTargetUptList))) / 20;
end
solution = optimizeCbModel(modelUpt);
if solution.stat~=1
    error('modelUpt model does not solve properly')
end

% GBA deletion 100%
gbaIdx = model.genes(contains('2629.', model.genes));
[GBA, ~, ~, ~] = deleteModelGenes(model, gbaIdx, 0);
    
%% In silico vs in vitro
objective = 'unWeightedTCBMfluxConc';
switch objective
    case 'unWeightedTCBMflux'
        % unWeightedTCBMflux
        model = modelUpt;
        model.osenseStr = 'min';
        model.cf = 0;
        model.cr = 0;
        model.g = 2;
        model.u0 = 0;
        model.f = 1;
        tcbmParam.method = 'fluxes';
        tcbmParam.printLevel = 1;
        tcbmParam.solver = 'mosek';
        [ctrlsol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
        
        % GBA deletion 100%
        [GBA, ~, ~, ~] = deleteModelGenes(model, gbaIdx, 0);
        [GBA_unWeightedTCBMflux, ~] = entropicFluxBalanceAnalysis(GBA, tcbmParam);
        GBAsol = GBA_unWeightedTCBMflux;
        
    case 'unWeightedTCBMfluxConc'
        % unWeightedTCBMfluxConc
        model = modelUpt;
        model.osenseStr = 'min';
        model.cf = 0;
        model.cr = 0;
        model.g = 2;
        model.u0 = 0;
        model.f = 1;
        tcbmParam.method = 'fluxConc';
        tcbmParam.printLevel = 1;
        tcbmParam.solver = 'mosek';
        [ctrlsol, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
        
        % GBA deletion 100%
        [GBA, ~, ~, ~] = deleteModelGenes(model, gbaIdx, 0);
        [GBA_unWeightedTCBMflux, ~] = entropicFluxBalanceAnalysis(GBA, tcbmParam);
        GBAsol = GBA_unWeightedTCBMflux;
end
predSol1 = ctrlsol;
predSol1.name = 'control';
predSol2 = GBAsol;
predSol2.name = 'GBA1 ko';




labelsL = strtrim(exoMet.name);
labelsR = exoMet.platform;
       
nlt=length(labelsL);
bool=true(nlt,1);
maxCharacters=0;
for i=1:nlt
    nCharacters=length(labelsL{i});
    if nCharacters>maxCharacters
        maxCharacters=nCharacters;
    end
    if isempty(labelsL{i})
        labelsL{i}='';
        bool(i)=0;
    else
        if ismember(labelsL{i},{'Nc20:4 ';'Ex_icit[e]';'Nc20:4 ';'EX_CE1557[e]';'EX_HC00900[e]';'Transport of (R)-3-Hydroxybutanoate  via H+ Symport';'EX_M03117[e]'})
            bool(i)=0;
        end
    end
    if isempty(labelsR{i})
        labelsR{i}='';
    end
end
labelsL = strrep(labelsL,'Exchange of ','');
labelsL = strrep(labelsL,'3, 4-dihydroxy','3,4-dihydroxy');
labelsL = strrep(labelsL,'N-Acetyl-Tyrosine','N-Acetyl-L-tyrosine');
labelsL = strtrim(labelsL);
exoMet.labelsL=labelsL;

labelsR = strrep(labelsR,'GC-MS','$GCMS$');
labelsR = strrep(labelsR,'LC-MS','$AccQtag$');
labelsR = strrep(labelsR,'dmpa_old','$DmPa^{13}C$');
labelsR = strrep(labelsR,'dmpa_new','$DmPa^{2}H$');
labelsR = strrep(labelsR,'dmpaN','$DmPa^{2}H$');
labelsR = strrep(labelsR,'dmpa','$DmPa^{13}C$');
labelsR = strrep(labelsR,'bzcl_old','$BzCl$');
labelsR = strrep(labelsR,'bzcl_new','$BzCl$');
labelsR = strrep(labelsR,'bzcl','$BzCl$');
exoMet.labelsR=labelsR;


switch plotParam.labelType
    case 'metabolitePlatform'
        for i=1:nlt
            numberOfCharacters = maxCharacters - length(labelsL{i}) + 1;
            numberOfCharacters = 1;
            labels{i}=[labelsL{i} ' ' repmat('\ ',1,numberOfCharacters) ' ' labelsR{i}];
        end
    case 'metabolite'
        labels = labelsL;
    case 'platform'
        labels = labelsR;
end
exoMet.labels=labels;

%% plotExperimentalvsPredictedExchange

% Measurements vs model
plotExperimentalvsPredictedExchange(model, exoMet, predSol1, plotParam)

% Model vs C1 inhibition
C1 = changeRxnBounds(model, 'NADH2_u10mi', 0, 'b');

C1 = modelUpt;
C1.osenseStr = 'min';
C1.cf = 0;
C1.cr = 0;
C1.g = 2;
C1.u0 = 0;
C1.f = 1;
tcbmParam.method = 'fluxConc';
tcbmParam.printLevel = 1;
tcbmParam.solver = 'mosek';
[C1sol, ~] = entropicFluxBalanceAnalysis(C1, tcbmParam);
predSol1 = ctrlsol;
predSol1.name = 'control';
predSol2 = C1sol;
predSol2.name = 'δv (prediction)';

% absolute change in flux
exoMetData = '~/work/sbgCloud/programExperimental/projects/PINK1SysBio/data/metabolomicsData/FinalData/summarisedData/combined_platforms_rates.mat';
load(exoMetData)
absoluteFluxChange = table;
absoluteFluxChange.rxnID = combinedData.rxnID;
absoluteFluxChange.rxnNames = strcat('Exchange of ',combinedData.metName);
absoluteFluxChange.mean = combinedData.mean_ctrl - combinedData.mean_rot;
absoluteFluxChange.SD = (combinedData.sd_ctrl + combinedData.sd_rot) / 2;
absoluteFluxChange.units(:) = {'μmol/gDW/hr'};
absoluteFluxChange.platform = combinedData.platform;
[~,rxnID] = ismember(exoMet.rxnID, absoluteFluxChange.rxnID);
absoluteFluxChange.labels(rxnID) = exoMet.labels;
absoluteFluxChange.changedSing = sign(combinedData.mean_ctrl) .* sign(combinedData.mean_rot) == -1;
toDelete = ~ismember(absoluteFluxChange.rxnID, model.rxns) | isnan(absoluteFluxChange.SD);
absoluteFluxChange(toDelete, :) = [];

plotExperimentalvsPredictedExchange(model, absoluteFluxChange, predSol2)
title('Comparison of measured and predicted changes in reaction rates, δ v, between control and complex I inhibition')

% Correlation
combinedData2 = combinedData;
toDelete = ~ismember(combinedData2.rxnID, model.rxns) | isnan(combinedData2.sd_rot);
combinedData2(toDelete, :) = [];

display(['Spearman correlation: ' num2str(corr(predSol2.v(findRxnIDs(model, combinedData2.rxnID)), combinedData2.mean_rot,'Type','Spearman'))])
