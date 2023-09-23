function transcriptomicData = transcriptomics2genes(file)
% Function for processed gene expression-based model generation
%
% USAGE:
%
%    [OmicsMapped] = findGenesInModelFromTrascriptomicAnalysis(model, file, varyThresh, nThresh)
%
% INPUTS:

%    file:          Path to the directory containing the transcriptomic
%                   analysis.
%
%
% OUTPUTS:
%    OmicsMapped:   Struct file with the genes indentified in the
%                   trascriptomic analysis in a COBRA model 
%
% .. Author: - Jennifer Modamio 08/06/2016

%% 2. Read and load data

% Read the data file
data = regexp(fileread(file), '\n', 'split')';

% If the last row is empty delete it
if ~isempty(data{end})
    dataLength = length(data);
else
    dataLength = length(data)-1;
end

%% 3. Convert all the samples expression into a matrix and the EntrezIDs in a string array

counter = 0;
for i = 1:dataLength
    dLine = strsplit(data{i}, '\t');
    entrezIDsString = strsplit(dLine{15}, ',');
    samples(i, 1:11) = cellfun(@str2double, {dLine{3:13}});
    entrezIDs{i} = entrezIDsString{1};
    % If there is two or more EntrezIDs repeat all the samples expression
    % at the end of the matrix
    if length(entrezIDsString)>1
        for j = 2:length(entrezIDsString)
            counter = counter + 1;
            samples(dataLength + counter, 1:11) = cellfun(@str2double, {dLine{3:13}});
            entrezIDs{dataLength+counter} = entrezIDsString{j};
        end
    end
end

uniqueEntrezIDs = unique(entrezIDs);

%% 4. Merge data
% To obtain one expression value per gene we are taking the max expression value
% for repeated genes and the mean of all samples (replicates) per gene

% Take the mean of the maximum values of the repeated gene expression values in the
% samples and the mean of non repeated

counter = 0;
for i = 1:length(uniqueEntrezIDs)
    entrezID = uniqueEntrezIDs{i};
    if ~isempty(entrezID) && ~isequal(entrezID, 'EntrezID') && ~isequal(entrezID, 'NA')
        expressionValues = samples(strcmp(entrezID, entrezIDs), :);
        [rows, ~]=size(expressionValues);
        if rows == 1
            expressionValuesMean = mean(expressionValues(expressionValues~=0));
            if expressionValuesMean ~= 0 && ~isnan(expressionValuesMean)
                counter = counter + 1;
                geneList{counter} = entrezID;
                expressionMean(counter) = expressionValuesMean;
            end
        else
            expressionValues = max(expressionValues);
            expressionValuesMean = mean(expressionValues(expressionValues~=0));
            if expressionValuesMean ~= 0 && ~isnan(expressionValuesMean)
                counter = counter + 1;
                geneList{counter} = entrezID;
                expressionMean(counter) = expressionValuesMean;
            end
        end
    end
end

transcriptomicData = table(geneList', expressionMean','VariableNames',{'genes', 'expVal'});