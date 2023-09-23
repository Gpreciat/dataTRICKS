function modelStat = modelData(model)
% Print the cobra model's statistics
%
% USAGE:
%
%    modelStat = modelData(model)
%
% INPUT:
%    model:     A generic COBRA model
%
%                     * .S - Stoichiometric matrix
%                     * .mets - Metabolite ID vector
%                     * .rxns - Reaction ID vector
%                     * .lb - Lower bound vector
%                     * .ub - Upper bound vector
%
% OUTPUT:
%    modelStat:  Table the the model statistics
%

if isfield(model, 'SIntRxnBool')
    nonInternalRxnIdxs =  find(~model.SIntRxnBool);
else
    nonInternalRxnIdxs = strncmp('EX_', model.rxns, 3) == 1;
    nonInternalRxnIdxs = nonInternalRxnIdxs + strncmp('DM_', model.rxns, 3) == 1;
    nonInternalRxnIdxs = find(nonInternalRxnIdxs + strncmp('sink_', model.rxns, 5) == 1);    
end

closedEx = [];
secEx = [];
uptEx = [];
revEx = [];
for i = 1:length(nonInternalRxnIdxs)
    if  model.lb(nonInternalRxnIdxs(i)) < 0 && model.ub(nonInternalRxnIdxs(i)) <= 0
        uptEx = [uptEx nonInternalRxnIdxs(i)];
    elseif model.lb(nonInternalRxnIdxs(i)) >= 0 && model.ub(nonInternalRxnIdxs(i)) > 0
        secEx = [secEx nonInternalRxnIdxs(i)];
    elseif model.lb(nonInternalRxnIdxs(i)) < 0 && model.ub(nonInternalRxnIdxs(i)) > 0
        revEx = [revEx nonInternalRxnIdxs(i)];
    else
        closedEx = [closedEx nonInternalRxnIdxs(i)];
    end
end

S = full(model.S);
[mlt, nlt] = size(S);

modelInfo = table([...
    mlt; ...
    numel(unique(regexprep(model.mets, '(\[\w\])', ''))); ...
    nlt; ...
    numel(setdiff(1:nlt, nonInternalRxnIdxs)); ...
    numel(model.rxns(nonInternalRxnIdxs)); ...
    numel(closedEx); ...
    numel(secEx); ...
    numel(uptEx); ...
    numel(revEx); ...
    numel(model.genes); ...
    rank(full(model.S)); ...
    numel(unique([model.subSystems{:}]'))], ...
    ...
    'VariableNames', ...
    {'Statistics'},...
    'RowNames',...
    {'Metabolites'; ...
    'Unique metabolites'; ...
    'Reactions'; ...
    'Internal reactions'; ...
    'External reactions'
    'External reactions closed'; ...
    'External reactions secretions'; ...
    'External reactions uptakes'; ...
    'External reactions reversibles'; ...
    'Genes'; ...
    'rank(S)'; ...
    'Sub systems'});
% 
% modelInfo.rxns = model.rxns;
% modelInfo.internalRxns = model.rxns;
% modelInfo.internalRxns(nonInternalRxnIdxs) = [];
% modelInfo.subSystems = unique([model.subSystems{:}]');
% modelInfo.exRxns = model.rxns(nonInternalRxnIdxs);
% modelInfo.closedEx = model.rxns(closedEx);
% modelInfo.secEx = model.rxns(secEx);
% modelInfo.uptEx = model.rxns(uptEx);
% modelInfo.revEx = model.rxns(revEx);
% modelInfo.mets = model.mets;
% modelInfo.uniqueMets = unique(regexprep(model.mets, '(\[\w\])', ''));
% modelInfo.genes = model.genes;
% modelInfo.sparsity = (numel(find(S)) *100)/numel(S);
