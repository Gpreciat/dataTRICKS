function [rxnTable, metTable] = commonIdsBetweenModels(modelBase, modelReference, idType, checkCompartments)
% Given two input models, it compares their metabolite and reaction identifiers
% and returns the set of identifiers that are present in both models. This
% function is useful for reconciling and combining different metabolic network
% models, as it allows users to find overlapping components in order to
% construct a unified and more comprehensive network.
%
% USAGE:
%
%    [rxnTable, metTable] = commonIdsBetweenModels(modelBase, modelReference, idType, checkCompartments)
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
    matchingRxnsFields = {'rxnkegg'; 'rxnbiocyc'; 'rxnmetanet'; 'rxnrhea'; ...
        'rxnbigg'};
elseif isequal(idType, 'both')
    matchingMetsFields = {'metkegg'; 'metchebi'; 'metinchi'; 'metsmiles'; 'methmdb'; ...
        'metpubchem'; 'metseed'; 'metbigg'; 'metbiocyc'; 'metNames'};
    matchingRxnsFields = {'rxnkegg'; 'rxnbiocyc'; 'rxnmetanet'; 'rxnrhea'; ...
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
%     modelNew = alphabetizeModel(modelNew);
%     modelNew.rxnFormulas = printRxnFormula(modelNew, 'printFlag', 0);
%     modelReference = alphabetizeModel(modelReference);
%     modelReference.rxnFormulas = printRxnFormula(modelReference, 'printFlag', 0);
    
    % Make a table with the information of all the metabolites in both
    % models
    nRows = length(modelNew.mets) + sum(~ismember(modelReference.mets, modelNew.mets));
    varTypes = {'string', 'string', 'string', 'string', 'string', 'double', ...
        'string', 'string', 'string', 'string', 'string', 'string', 'string',...
        'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string',...
        'string', 'string', 'string'};
    varNames = {'mets', 'status', 'comps', 'metNames', 'metFormulas', ...
        'metCharges', 'metHMDBIDbase', 'metHMDBIDref', 'metInChIStringbase', ...
        'metInChIStringref', 'metKEGGIDbase', 'metKEGGIDref', 'metChEBIIDbase', ...
        'metChEBIIDref', 'metMetaNetXIDbase', 'metMetaNetXIDref', 'metSEEDIDbase', ...
        'metSEEDIDref', 'metBiGGIDbase', 'metBiGGIDref', 'metBioCycIDbase', ...
        'metBioCycIDref', 'metSBOTermsbase', 'metSBOTermsref', 'inchikeybase', ...
        'inchikeyref'};
    metTable = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
        'VariableNames', varNames);
    
    % modelNew
    for i = 1:length(modelNew.mets)
        metTable.mets(i) = modelNew.mets(i);
        if ismember(modelNew.mets(i), modelReference.mets)
            metTable.status(i) = 'both';
        else
            metTable.status(i) = 'onlyBase';
        end
        metTable.comps(i) = regexprep(modelNew.mets(i), '.*\[(.*?)\]', '$1');
        if any(contains(modelBaseFields, 'metnames', IgnoreCase=true))
            metTable.metNames(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metnames', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metformulas', IgnoreCase=true))
            metTable.metFormulas(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metformulas', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metcharge', IgnoreCase=true))
            metTable.metCharges(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metcharge', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'methmdb', IgnoreCase=true))
            metTable.metHMDBIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'methmdb', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'inchis', IgnoreCase=true))
            metTable.metInChIStringbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'inchis', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metkegg', IgnoreCase=true))
            metTable.metKEGGIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metkegg', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metchebi', IgnoreCase=true))
            metTable.metChEBIIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metchebi', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metmetanet', IgnoreCase=true))
            metTable.metMetaNetXIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metmetanet', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metseed', IgnoreCase=true))
            metTable.metSEEDIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metseed', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metbigg', IgnoreCase=true))
            metTable.metBiGGIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metbigg', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metbiocyc', IgnoreCase=true))
            metTable.metBioCycIDbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metbiocyc', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'metsbo', IgnoreCase=true))
            metTable.metSBOTermsbase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'metsbo', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'inchikey', IgnoreCase=true))
            metTable.inchikeybase(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'inchikey', IgnoreCase=true)})(i);
        end
    end
    
    % modelNew
    idx = 0;
    for i = 1:length(modelReference.mets)
        idBool = ismember(metTable.mets, modelReference.mets(i));
        if ~any(idBool)
            idx = idx + 1;
            idBool(length(modelNew.mets) + idx) = true;
            metTable.mets(idBool) = modelReference.mets(i);
            metTable.status(idBool) = 'onlyRef';
            metTable.comps(idBool) = regexprep(modelReference.mets(i), '.*\[(.*?)\]', '$1');
            metTable.metNames(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metnames', IgnoreCase=true)})(i);
            metTable.metFormulas(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metformulas', IgnoreCase=true)})(i);
            if any(contains(modelReferenceFields, 'metcharge', IgnoreCase=true))
                metTable.metCharges(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metcharge', IgnoreCase=true)})(i);
            end
        end
        if any(contains(modelReferenceFields, 'methmdb', IgnoreCase=true))
            metTable.metHMDBIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'methmdb', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'inchis', IgnoreCase=true))
            metTable.metInChIStringref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'inchis', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metkegg', IgnoreCase=true))
            metTable.metKEGGIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metkegg', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metchebi', IgnoreCase=true))
            metTable.metChEBIIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metchebi', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metmetanet', IgnoreCase=true))
            metTable.metMetaNetXIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metmetanet', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metseed', IgnoreCase=true))
            metTable.metSEEDIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metseed', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metbigg', IgnoreCase=true))
            metTable.metBiGGIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metbigg', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metbiocyc', IgnoreCase=true))
            metTable.metBioCycIDref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metbiocyc', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'metsbo', IgnoreCase=true))
            metTable.metSBOTermsref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'metsbo', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'inchikey', IgnoreCase=true))
            metTable.inchikeyref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'inchikey', IgnoreCase=true)})(i);
        end
    end
    
end

if exist('matchingRxnsFields', 'var')
    
    % Check wich fields are shared between both models
    matchingRxnsFields = matchFields(modelBaseFields, modelReferenceFields, matchingRxnsFields);
    
    % Identifies from each of the matchingFields, the same id
    rxnsInCommon = findSameIds(modelBase, modelReference, matchingRxnsFields, 'rxns');
    
    % Rename the metabolites with the same identifier
    for i = 1:length(rxnsInCommon)
        if ~isempty(rxnsInCommon{i, 2})
            modelNew.rxns(rxnsInCommon{i, 2}) = modelReference.rxns(rxnsInCommon{i, 1});
        end
    end
    modelNew = alphabetizeModel(modelNew);
    modelNew.rxnFormulas = printRxnFormula(modelNew, 'printFlag', 0);
    modelReference = alphabetizeModel(modelReference);
    modelReference.rxnFormulas = printRxnFormula(modelReference, 'printFlag', 0);
    
    % Make a table with the information of all the reactions in both
    % models
    nRows = length(modelNew.rxns) + sum(~ismember(modelReference.rxns, modelNew.rxns));
    varTypes = {'string', 'string', 'double', 'double', 'double', 'double', ...
        'string', 'string', 'string', 'string', 'string', 'string', 'string', ...
        'string', 'string', 'string', 'string', 'string', 'string', 'string', ...
        'string', 'string', 'string', 'string', 'string', 'string', 'string', ...
        'string'};
    varNames = {'rxns', 'status', 'lb_base', 'lb_ref', 'ub_base','ub_ref',...
        'rxnFormulas_base','rxnFormulas_ref', 'rxnNames_base', 'rxnNames_ref', ...
        'rxnECNumbers_base', 'rxnECNumbers_ref', 'grRules_base', 'grRules_ref', ...
        'rules_base', 'rules_ref', 'rxnKEGGID_base', 'rxnKEGGID_ref', 'rxnMetaNetXID_base',...
        'rxnMetaNetXID_ref', 'rxnBioCycID_base', 'rxnBioCycID_ref','rxnRheaID_base',...
        'rxnRheaID_ref', 'rxnBiGGID_base', 'rxnBiGGID_ref', 'rxnSBOTerms_base',...
        'rxnSBOTerms_ref'};
    rxnTable = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes,...
        'VariableNames', varNames);
    
    matchingRxnsFields = {'ec'; 'rxnkegg'; 'rxnbiocyc'; 'rxnmetanet'; 'rxnrhea'; ...
        'rxnbigg'};
    
    % modelNew
    for i = 1:length(modelNew.rxns)
        rxnTable.rxns(i) = modelNew.rxns(i);
        if ismember(modelNew.rxns(i), modelReference.rxns)
            rxnTable.status(i) = 'both';
        else
            rxnTable.status(i) = 'onlyBase';
        end
        rxnTable.lb_base(i) = modelNew.lb(i);
        rxnTable.ub_base(i) = modelNew.ub(i);
        rxnTable.rxnFormulas_base(i) = modelNew.rxnFormulas(i);
        rxnTable.rxnNames_base(i) = modelNew.rxnNames(i);
        if any(contains(modelBaseFields, 'ec', IgnoreCase=true))
            rxnTable.rxnECNumbers_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'ec', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'grRules', IgnoreCase=true))
            rxnTable.grRules_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'grRules', IgnoreCase=true)})(i);
        end
        if ismember(modelBaseFields, 'rules')
            rxnTable.rules_base(i) = modelNew.(modelBaseFields{ismember(modelBaseFields, 'rules')})(i);
        end
        if any(contains(modelBaseFields, 'rxnkegg', IgnoreCase=true))
            rxnTable.rxnKEGGID_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'rxnkegg', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'rxnmetanet', IgnoreCase=true))
            rxnTable.rxnMetaNetXID_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'rxnmetanet', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'rxnbiocyc', IgnoreCase=true))
            rxnTable.rxnBioCycID_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'rxnbiocyc', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'rxnrhea', IgnoreCase=true))
            rxnTable.rxnRheaID_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'rxnrhea', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'rxnbigg', IgnoreCase=true))
            rxnTable.rxnBiGGID_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'rxnbigg', IgnoreCase=true)})(i);
        end
        if any(contains(modelBaseFields, 'rxnsbo', IgnoreCase=true))
            rxnTable.rxnSBOTerms_base(i) = modelNew.(modelBaseFields{contains(modelBaseFields, 'rxnsbo', IgnoreCase=true)})(i);
        end
    end
    
    % modelNew
    idx = 0;
    for i = 1:length(modelReference.rxns)
        idBool = ismember(rxnTable.rxns, modelReference.rxns(i));
        if ~any(idBool)
            idx = idx + 1;
            idBool(length(modelNew.rxns) + idx) = true;
            rxnTable.rxns(idBool) = modelReference.rxns(i);
            rxnTable.status(idBool) = 'onlyRef';           
        end
        rxnTable.lb_ref(idBool) = modelReference.lb(i);
        rxnTable.ub_ref(idBool) = modelReference.ub(i);
        rxnTable.rxnFormulas_ref(idBool) = modelReference.rxnFormulas(i);
        rxnTable.rxnNames_ref(idBool) = modelReference.rxnNames(i);
        if any(contains(modelReferenceFields, 'ec', IgnoreCase=true))
            rxnTable.rxnECNumbers_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'ec', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'grRules', IgnoreCase=true))
            rxnTable.grRules_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'grRules', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rules', IgnoreCase=true))
            rxnTable.rules_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rules', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rxnkegg', IgnoreCase=true))
            rxnTable.rxnKEGGID_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rxnkegg', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rxnmetanet', IgnoreCase=true))
            rxnTable.rxnMetaNetXID_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rxnmetanet', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rxnbiocyc', IgnoreCase=true))
            rxnTable.rxnBioCycID_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rxnbiocyc', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rxnrhea', IgnoreCase=true))
            rxnTable.rxnRheaID_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rxnrhea', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rxnbigg', IgnoreCase=true))
            rxnTable.rxnBiGGID_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rxnbigg', IgnoreCase=true)})(i);
        end
        if any(contains(modelReferenceFields, 'rxnsbo', IgnoreCase=true))
            rxnTable.rxnSBOTerms_ref(idBool) = modelReference.(modelReferenceFields{contains(modelReferenceFields, 'rxnsbo', IgnoreCase=true)})(i);
        end
        
    end
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