%identify the steps that remove active genes from the model

%generic model
inputFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programExperimental' ...
    filesep 'projects' filesep 'xomics' filesep 'data' filesep 'Recon3D_301'];
genericModelName = 'Recon3DModel_301_xomics_input.mat';
load([inputFolder filesep genericModelName])
genericGenesStr = model.genes;
genericGenesStr = regexprep(genericGenesStr, '\.\d', '');
genericGenesStr = cellfun(@str2num,genericGenesStr,'UniformOutput',false);
genericGenes=NaN*ones(length(genericGenesStr),1);
for i=1:length(genericGenesStr)
    genericGenes(i,1) = genericGenesStr{i};
end

% specificData
dataFolder = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'data' filesep 'xomics' filesep];
bibliomicData = 'specificData.xlsx';
specificData = preprocessingOmicsModel([dataFolder bibliomicData], 1, 1);
%specificData.exoMet = readtable([dataFolder 'exoMet']);
%specificData.transcriptomicData = readtable([dataFolder 'transcriptomicData']);
%specificData.transcriptomicData.genes = string(specificData.transcriptomicData.genes);
activeGenes = specificData.activeGenes;

cd ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/iDopaNeuro1
debugFiles = dir;
debugFiles = {debugFiles.name};
debugFiles = debugFiles(contains(debugFiles, 'debug'))';
debugFilesNum = strtok(debugFiles,'.');
debugFilesNum = cellfun(@str2num,debugFilesNum);
[~,xi]=sortrows(debugFilesNum,'ascend');
debugFiles = debugFiles(xi);

fprintf('%13s %s      %s\n','#Active_genes','         ','Stage')
fprintf('%13u %13u  %s\n',length(activeGenes),NaN,'Active gene list')

fprintf('%13s %s  %s\n','#Active_genes','#Model_genes','Stage')
bool = ismember(activeGenes,genericGenes);
fprintf('%13u %13u %s\n',nnz(bool),length(isfinite(genericGenes)),'Generic model')


for j=1:length(debugFiles)
    %load the model at each stage
    load(debugFiles{j});
    
    genericGenesStr = model.genes;
    genericGenesStr = regexprep(genericGenesStr, '\.\d', '');
    genericGenesStr = cellfun(@str2num,genericGenesStr,'UniformOutput',false);
    genericGenes=NaN*ones(length(genericGenesStr),1);
    for i=1:length(genericGenesStr)
        genericGenes(i,1) = genericGenesStr{i};
    end
    
    bool = ismember(activeGenes,genericGenes);
    fprintf('%13u %13u %s\n',nnz(bool),length(isfinite(genericGenes)), debugFiles{j})
end

