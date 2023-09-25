function modelNew = commonIdsBetweenModels(modelBase, modelReference, idType, checkCompartments)
% Given two input models, it compares their metabolite and reaction identifiers
% and returns the set of identifiers that are present in both models. This
% function is useful for reconciling and combining different metabolic network
% models, as it allows users to find overlapping components in order to
% construct a unified and more comprehensive network.
%
% USAGE:
%
%    modelNew = commonIdsBetweenModels(modelBase, modelReference, idType, checkCompartments)
%
% INPUT:
%   modelBase:  Base COBRA model with identifiers.
%   modelReference:  Is a reference COBRA model used to find common
%               identifiers in "modelBase".
%
% OPTIONAL INPUTS:
%   idType:     String indicating if the ids to correlate come from
%               metabolites ('mets'), reactions ('rxns') or both ('both')
%               (default: 'both').
%   checkCompartments: Boolean value to indicate whether the IDs will be checked
%               for the compartments as well (default: false).
%
% OUTPUT:
%   modelNew:   A new COBRA model with correlated identifiers from both
%               models.
%
% 2023 German Preciat

if (nargin < 3)
    idType = 'mets';
end
if (nargin < 4)
    checkCompartments = false;
end

% Create a new model in which ids will be modified
modelNew = modelBase;

% Identify the fields to compare
if isequal(idType, 'mets')
    matchingMetsFields = {'metkegg'; 'metchebi'; 'metinchi'; 'metsmiles'; 'methmdb'; ...
        'metpubchem'; 'metseed'; 'metbigg'; 'metbiocyc'; 'metNames'};
elseif isequal(idType, 'rxns')
    matchingRxnsFields = {'ec'; 'rxnkegg'; 'rxnbiocyc'; 'rxnmetanet'; 'rxnrhea'; ...
        'rxnbigg'};
elseif isequal(idType, 'both')
    matchingMetsFields = {'metkegg'; 'metchebi'; 'metinchi'; 'metsmiles'; 'methmdb'; ...
        'metpubchem'; 'metseed'; 'metbigg'; 'metbiocyc'; 'metNames'};
    matchingRxnsFields = {'ec'; 'rxnkegg'; 'rxnbiocyc'; 'rxnmetanet'; 'rxnrhea'; ...
        'rxnbigg'};
else
    error(['The idType ''' idType ''' is not recognised'])
end

% Extract the field name for each model
modelBaseFields = fieldnames(modelBase);
modelReferenceFields = fieldnames(modelReference);

if exist('matchingMetsFields', 'var')
    
    % Check wich fields are shared between both models
    matchingMetsFields = matchFields(modelBaseFields, modelReferenceFields, matchingMetsFields);
    
    % Identifies from each of the matchingFields, the same id
    metsInCommon = findSameIds(modelBase, modelReference, matchingMetsFields, 'mets');
    
    % Rename the metabolites with the same identifier
    for i = 1:length(metsInCommon)
        if ~isempty(metsInCommon{i, 2})
            modelNew.mets(metsInCommon{i, 2}) = strcat(...
                regexprep(modelReference.mets(metsInCommon{i, 1}), '(\[\w\])', ''), ...
                '[', ...
                regexprep(modelNew.mets(metsInCommon{i, 2}), '.*\[(.*?)\]', '$1'),...
                ']');
        end
    end
    modelNew = alphabetizeModel(modelNew);
    modelNew.rxnFormulas = printRxnFormula(modelNew, 'printFlag', 0);
        
end

% Find matching formulas
modelNewFormulas = regexprep(modelNew.rxnFormulas, '(\[\w\])', '');
modelReference = alphabetizeModel(modelReference);
modelReference.rxnFormulas = printRxnFormula(modelReference, 'printFlag', 0);
modelReferenceFormulas = regexprep(modelReference.rxnFormulas, '(\[\w\])', '');
sum(ismember(modelReferenceFormulas, modelNewFormulas))


if exist('matchingRxnsFields', 'var')
    
    % Check wich fields are shared between both models
    matchingRxnsFields = matchFields(modelBaseFields, modelReferenceFields, matchingRxnsFields);
    
    % Identifies which fields in both models have the same identifiers
    rxnsInCommon = findSameIds(modelBase, modelReference, matchingRxnsFields, 'rxns');
    
    
    
end

end

function matchingFields = matchFields(modelBaseFields, modelReferenceFields, listOfFields)
% Check wich fields are shared between both models

% Identify common IDs
fieldsBool = false(size(listOfFields));
for i = 1:length(listOfFields)
    if any(contains(modelBaseFields, listOfFields{i}, IgnoreCase=true)) && ...
            any(contains(modelReferenceFields, listOfFields{i}, IgnoreCase=true))
        fieldsBool(i) = true;
    end
end
matchingFields = listOfFields(fieldsBool);

end

function commonIds = findSameIds(modelBase, modelReference, matchingFields, idType)
% Identifies from each of the matchingFields, the same id

modelBaseFields = fieldnames(modelBase);
modelReferenceFields = fieldnames(modelReference);

commonIds = cell(length(modelReference.(idType)), 2);
commonIds(:, 1) = num2cell(1:length(modelReference.(idType)));
for i = 1:length(matchingFields)
    referenceList = modelReference.(modelReferenceFields{find(contains(modelReferenceFields, matchingFields{i}, IgnoreCase=true))});
    baseList = modelBase.(modelBaseFields{find(contains(modelBaseFields, matchingFields{i}, IgnoreCase=true))});
    for j = 1:length(referenceList)
        ids = split(referenceList{j}, '; ');
        for k = 1:length(ids)
            if any(contains(baseList, ids{k})) && ~isempty(ids{k})
                commonIds{j, 2} = unique([commonIds{j, 2}; find(contains(baseList, ids{k}))]);
            end
        end
    end
end

end