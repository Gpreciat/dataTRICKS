% Comparison of the fraction of reactions in (red) and out (blue) of the minimal
% flux vector in each metabolic subsystem of the iDopaNeuro1 model, obtained by
% minimising the function Ψ(v):=||v||_{0} subject to v∈Ω,
% to predict the minimum number of reactions that are required to be active to satisfy
% dopaminergic neuron specific constraints on the steady-state flux space Ω.
%
% The graph is saved in:
% ~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/codeResults/iDN1

clear

% define directories
pathSave = ['~' filesep 'work' filesep 'sbgCloud' filesep 'programReconstruction' ...
    filesep 'projects' filesep 'exoMetDN' filesep 'results' filesep 'codeResults' filesep 'iDN1'];
[~, ~] = changeCobraSolver('gurobi', 'all', 0);

models = {'iDopaNeuroCT'; 'iDopaNeuroC'};
for j = 1:length(models)
    
    % Load model
    load([pathSave filesep models{j} filesep models{j} '.mat'])
    eval(['model = ' models{j} ';'])
    model.subSystems(cellfun(@iscell, model.subSystems)) = [model.subSystems{cellfun(@iscell, model.subSystems)}];
    
    mediumRxns = model.XomicsToModelSpecificData.mediaData.rxnID;
    coreRxnAbbr = model.rxns(model.activeInactiveRxn == 1);
    
    % sparseFBA parameters
    changeCobraSolver('gurobi', 'LP');
    osenseStr = 'none';
    checkMinimalSet = true;
    checkEssentialSet = true;
    zeroNormApprox  = 'all';
    TolSparse = getCobraSolverParams('LP', 'feasTol');
    
    modelCoupled = model;
    iDopaNeuro2 = model;
    
    if ~isfield(modelCoupled, 'C')
        model.C = [];
    else
        nd = iDopaNeuro2.d;
        for i = 1:length(nd)
            ncsense = length(iDopaNeuro2.csense);
            iDopaNeuro2.csense((ncsense + 1), 1) = 'G';
        end
    end
    
    % FBA gives the minimum reactionBool of active reactions
    [v, rxnbool, essentialBool]  = sparseFBA(iDopaNeuro2, osenseStr,...
        checkMinimalSet,checkEssentialSet, zeroNormApprox);
    % ---Non-convex approximation---
    % 393 of these are heuristically minimal rxns.
    % 392 of these are essential rxns.
    
    % find heuristically minimal
    essential = model.rxns(find(essentialBool));
    heuristically = unique([model.rxns(find(rxnbool)); essential]);
    
    % find the minimum media composition
    rx_minMedia = intersect(essential, mediumRxns);
    %rx_heurMedia = intersect(heuristically, mediumRxns)
    %extra= setdiff(rx_heurMedia,rx_minMedia)
    summary_minMedia= [model.rxns(findRxnIDs(model,...
        rx_minMedia)) num2cell(model.lb(findRxnIDs(model,...
        rx_minMedia))) num2cell(model.ub(findRxnIDs(model,...
        rx_minMedia)))];
    %summary_minMedia2= [model_hNESC2DN.rxns(findRxnIDs(model_hNESC2DN,rx_heurMedia)) num2cell(model_hNESC2DN.lb(findRxnIDs(model_hNESC2DN,rx_heurMedia))) num2cell(model_hNESC2DN.ub(findRxnIDs(model_hNESC2DN,rx_heurMedia)))]
    summaryEssentialRxnBool = [model.rxns(findRxnIDs(model, essential)) ...
        num2cell(model.lb(findRxnIDs(model, essential))) ...
        num2cell(model.ub(findRxnIDs(model, essential))) ...
        model.subSystems(findRxnIDs(model, essential))];
    
    % excluded media metabolites
    id = findRxnIDs(model,setdiff(mediumRxns,heuristically));
    id = id(find(id));
    summary_excludedMedia= [model.rxns(id) num2cell(model.lb(id)) num2cell(model.ub(id))];
    
    % find how many from the CoreRxnlist
    heuristic_CoreRxns = intersect(findRxnIDs(model, coreRxnAbbr), findRxnIDs(model, heuristically)); %176 reactions from the essentialRxnsBool are in commom with CoreRxnList
    heuristic_CoreRxns_summary = [model.rxns(heuristic_CoreRxns) ...
        printRxnFormula(model, 'rxnAbbrList', model.rxns(heuristic_CoreRxns), 'printFlag', false) ...
        model.subSystems(heuristic_CoreRxns)];
    
    %% create inputs to visualise in the metabolic map
    % load map ReconMap3
    % % % %
    % % % % mapPath = '~/work/sbgCloud/programReconstruction/projects/metabolicCartography/results/2017_ReconMap/VMHmaps/';
    % % % % mapfile = 'ReconMap3.0_default_8225_newRxnMetaID.mat'
    % % % % load([mapPath mapfile])
    % % % % % default color of the map
    % % % % map_default = defaultLookMap(map);
    % % % % % change rxns to light grey
    % % % % map_grey = changeRxnColorAndWidth(map_default,map_default.rxnName,'LIGHTGRAY',2);
    % % % %
    % % % % % color rxns
    % % % % map1 = changeRxnColorAndWidth(map_grey,iDopaNeuro.rxns,'DARKRED',18); % all map
    % % % % %map2 = changeRxnColorAndWidth(map1,heuristically,'DARKSALMON',18) % all map
    % % % % map3 = changeRxnColorAndWidth(map1,essential,'MEDIUMSLATEBLUE',18); % rxnsSparseFBA
    
    % transfrom map to xml file to be open in cell designer
    % transformMap2XML(xml, map3, 'Recon3map_sparseFBA.xml')
    
    %% generate a plot (reactions per subsystems)
    
    subsystems = unique(model.subSystems);%obtain unique subsystems
    index = find(~cellfun(@isempty, subsystems)); %avoid error for missing subsystems
    subsystems = subsystems(index);
    
    % look for reactions per subsystem in model (except sparseFBA)
    rxnsNo = findRxnIDs(model, setdiff(model.rxns,essential));
    subNo = [model.rxns(rxnsNo) model.subSystems(rxnsNo)];
    
    % summaryEssentialRxnBool
    essentialSub = unique([summaryEssentialRxnBool{:,4}])';
    
    % length subsystems(find reactions in three lists)
    tableSubsystems = [];
    for i = 1:length(subsystems)
        
        tableSubsystems(i, 1) = sum(contains([model.subSystems{rxnsNo}], subsystems{i}));
        tableSubsystems(i, 2) = sum(contains([summaryEssentialRxnBool{:,4}], subsystems{i}));
        tableSubsystems(i, 3) = sum(tableSubsystems(i,:));
        
    end
    % obtain percentages for representation
    for i = 1:length(subsystems)
        percentage1(i) = (tableSubsystems(i,1)*100)/tableSubsystems(i,3);
        percentage2(i) = (tableSubsystems(i,2)*100)/tableSubsystems(i,3);
    end
    
    [~, idx] = sort(percentage1);
    
    %plot  (vertical)
    figure()
    subsystem_summary = [subsystems(idx) num2cell(percentage1(idx))' num2cell(percentage2(idx))']; %delete subsystmes not in consistent models
    
    b = barh(cell2mat(subsystem_summary(:,[2 3])),'stacked','DisplayName','percentage','LineStyle','none' );
    set(b(1),'DisplayName','TotalModel');
    set(b(2),'DisplayName','SparseFBA');
    set(gca,'YTick',1:length(subsystem_summary),'YTickLabel',subsystem_summary(:,1));
    title('Sparse reactions per subsystem', 'FontSize', 16)
    
%     saveas(gcf,[pathSave filesep models{j} filesep 'sparseFBA.fig']);
    
    clearvars -except pathSave j models
end