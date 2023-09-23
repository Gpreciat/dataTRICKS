% Generates a venn diagram that compares teh core reactions between the
% iDopaNeuroCT and iDopaNeuroC models.
%
% The graph is saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1

%% Define directories

clear

% Results dir
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];

% specificData dir
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
bibliomicData = 'bibliomicData.xlsx';

%% Core reactions

% Read activeReactions
specificData = preprocessingOmicsModel([dataFolder bibliomicData], 1, 1);
% coreRxnAbbr
coreRxnAbbr = {};
coreRxnAbbr = [coreRxnAbbr; specificData.activeReactions];
coreRxnAbbr = [coreRxnAbbr; specificData.rxns2constrain.rxnID];
for i = 1:length(specificData.coupledRxns.coupledRxnsList)
    coreRxnAbbr = [coreRxnAbbr; split(specificData.coupledRxns.coupledRxnsList{i}, ', ')];
end
coreRxnAbbr = [coreRxnAbbr; specificData.mediaData.rxnID];
coreRxnAbbr = unique(coreRxnAbbr);

%% Venn diagram

% Load models
load([pathSave filesep 'iDopaNeuroCT' filesep 'iDopaNeuroCT.mat'])
load([pathSave filesep 'iDopaNeuroC' filesep 'iDopaNeuroC.mat'])

iDopaNeuroCT_rxns = iDopaNeuroCT.rxns;
iDopaNeuroC_rxns = iDopaNeuroC.rxns;

% Generate diagram
set7 = numel(intersect(iDopaNeuroC_rxns, intersect(iDopaNeuroCT_rxns, coreRxnAbbr)));
vennX([numel(setdiff(coreRxnAbbr, [iDopaNeuroCT_rxns; iDopaNeuroC_rxns])) ...
    numel(intersect(coreRxnAbbr, iDopaNeuroCT_rxns))  - set7 ...
    numel(setdiff(iDopaNeuroCT_rxns, [coreRxnAbbr; iDopaNeuroC_rxns]))  ...
    numel(intersect(iDopaNeuroC_rxns, iDopaNeuroCT_rxns))  - set7 ...
    numel(setdiff(iDopaNeuroC_rxns, [iDopaNeuroCT_rxns; coreRxnAbbr]))  ...
    numel(intersect(coreRxnAbbr, iDopaNeuroC_rxns))  - set7 ...
    set7], .01);
title({'Core reactions between', 'iDopaNeuro models'}, 'FontSize', 16)

% Save figure
savefig([pathSave filesep models{j}  filesep 'vennMetabolomics'])
saveas(gcf,[pathSave filesep models{j}  filesep 'vennMetabolomics'],'png')
saveas(gcf,[pathSave filesep models{j}  filesep 'vennMetabolomics'],'eps')


