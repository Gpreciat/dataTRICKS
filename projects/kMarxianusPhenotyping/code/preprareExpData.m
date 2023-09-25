% This script is designed to process differential expression data for the 
% organism Kluyveromyces marxianus from a file and format the output for 
% compatibility with the COBRA Toolbox.
%
% Procedure:
% 
% Load Data: The script starts by reading the file 'diffExpr.txt'  
% containing differential expression data. The file is expected to have 
% columns such as "Transcript name", "Gene name", "logFC", "logCPM", "F", 
% "P-Value" and "FDR".
% 
% Data Transformation: The script performs necessary transformations, if 
% any, on the data. This may involve calculations or conversions to ensure 
% compatibility with COBRA Toolbox's requirements.
% 
% Formatting for COBRA: The processed data is then formatted into a 
% suitable structure that adheres to COBRA Toolbox's standards. This 
% includes organizing the data in a way that can be easily utilized for 
% subsequent analyses using COBRA Toolbox.
% 
% Output Generation: The formatted data is saved to an output file, 
% ensuring it is in a format that can be readily read and imported by 
% COBRA Toolbox. This output file becomes a bridge between the differential 
% expression analysis and the subsequent modeling steps.

%% Load Data

clear

% Define directories
filePath = regexprep(matlab.desktop.editor.getActiveFilename, ['code' filesep 'preprareExpData.m'], '');
dataDir = [filePath 'data' filesep];

% Read differential expression 
diffExp = readtable([dataDir 'diffExpr.txt']);

%% Formatting for COBRA

% Modify gene names
diffExp.GeneName = regexprep(diffExp.GeneName, '_(\d+)$', '');

% kmGEMv1 model
load([dataDir 'km.mat'])
kmGEMv1 = model;

% Change gene name for gene ID
for i = 1:length(diffExp.GeneName)
    geneIdBool = ismember(kmGEMv1.geneNames, diffExp.GeneName{i});
    if any(geneIdBool)
        diffExp.GeneName{i} = kmGEMv1.genes{geneIdBool};
    end
end

%% Data Transformation

kmExpData = table;
kmExpData.genes = diffExp.GeneName;
kmExpData.gene = diffExp.GeneName;
kmExpData.expVal = diffExp.logFC;
kmExpData.value = diffExp.logFC;

%% Output Generation

writetable(kmExpData, [dataDir 'kmExpData.txt'], 'QuoteStrings', false);