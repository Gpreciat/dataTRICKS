
% Load cornData
xlsDir = '~/work/sbgCloud/programReconstruction/projects/exoMetDN/data/omics/20210118_validationData_K7_final_ordered_umol_gDW_h.csv';
cornData = readtable(xlsDir);
cornRxns = cornData.exRxns;

% Load Can data
conc2Fluxes
canRxns = v(1, 2:end);
cornRxns(cellfun(@isempty, cornRxns)) = [];

% Check rxns
rxnsToCheck = intersect(canRxns, cornRxns);
rxnsToChange = setdiff(canRxns, cornRxns);

% Get the outlier
for i = 1:length(rxnsToCheck)
    idxCan = 1 + strmatch(rxnsToCheck{i}, v(1, 2:end), 'exact');
    idxCorn = strmatch(rxnsToCheck{i}, cornData.exRxns);
    can2cornData(i) = cornData.mean(idxCorn) / v{11, idxCan(1)};
end
can2CornOutlier = isoutlier(can2cornData);

figure
plot(can2cornData,[1:length(can2cornData)],'.')
yticks([1:length(can2cornData)])
set(gca,'TickLabelInterpreter','none')
labels=rxnsToCheck
yticklabels(labels)

factorScale = mean(can2cornData(~can2CornOutlier));

% factorScale = 0.479542184535408
%this approximately scales cans data to the scale of cornelius' data

for i = 1:length(rxnsToChange)
    idxCan = 1 + strmatch(rxnsToChange{i}, v(1, 2:end), 'exact');
    newValues(i) = v{11, idxCan(1)} * factorScale;
end
