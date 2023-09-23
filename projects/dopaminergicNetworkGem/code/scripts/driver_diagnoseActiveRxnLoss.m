%identify the steps that remove active genes from the model

%generic model
inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
    filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input.mat';
load([inputFolder filesep genericModelName])
genericRxns = model.rxns;

% specificData
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
bibliomicData = 'specificData.xlsx';
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

% Read dir
cd ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/iDopaNeuro1
debugFiles = dir;
debugFiles = {debugFiles.name};
debugFiles = debugFiles(contains(debugFiles, 'debug'))';
debugFilesNum = strtok(debugFiles,'.');
debugFilesNum = cellfun(@str2num,debugFilesNum);
[~,xi]=sortrows(debugFilesNum,'ascend');
debugFiles = debugFiles(xi);

fprintf('%13s %s      %s\n','#Active_rxns','         ','Stage')
fprintf('%13u %13u  %s\n',length(coreRxnAbbr),NaN,'Active gene list')

fprintf('%13s %s  %s\n','#Active_rxns','#Model_rxns','Stage')
bool = ismember(coreRxnAbbr,genericRxns);
fprintf('%13u %13u %s\n',nnz(bool),length(genericRxns),'Generic model')


for j=1:length(debugFiles)
    %load the model at each stage
    load(debugFiles{j});
    
    genericRxns = model.rxns;
    
    bool = ismember(coreRxnAbbr,genericRxns);
    fprintf('%13u %13u %s\n',nnz(bool),length(genericRxns), debugFiles{j})
end

