function [comparisonData] = quanQualAcc(model, predictedFlux, validationData, param)

%param.boundPrecisionLimit

if nargin < 4 || isempty(param)
    param = struct;
end

% Add default parameters
if ~isfield(param, 'boundPrecisionLimit')
    feasTol = getCobraSolverParams('LP', 'feasTol');
    param.boundPrecisionLimit = feasTol * 10;
end

if iscell(validationData.mean)
    validationData.mean = str2double(validationData.mean);
end
if iscell(validationData.SD)
    validationData.SD = str2double(validationData.SD);
end

% Experimental mean (experimentalMean < boundPrecisionLimit are considered as zero)
experimentalMean = validationData.mean;

% Identify secretions and uptakes present in the model
exoMetRxnsIdx = ismember(validationData.rxnID, model.rxns);

% Secretions
rxnTargetSecList = validationData.rxnID(experimentalMean > 0 & exoMetRxnsIdx);
rxnTargetSec_mean = validationData.mean(experimentalMean > 0 & exoMetRxnsIdx);
rxnTargetSec_SD = validationData.SD(experimentalMean > 0 & exoMetRxnsIdx);

% Uptakes
rxnTargetUptList = validationData.rxnID(experimentalMean < 0 & exoMetRxnsIdx);
rxnTargetUpt_mean = validationData.mean(experimentalMean < 0 & exoMetRxnsIdx);
rxnTargetUpt_SD = validationData.SD(experimentalMean < 0 & exoMetRxnsIdx);

% Prepare the table with the accuracy data
nRows = length(rxnTargetSecList) + length(rxnTargetUptList);
varTypes = {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
varNames = {'rxnID', 'mean', 'SD', 'target', 'v', 'predict', 'agree', 'dv'};
fullReport = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes, 'VariableNames', varNames);

% Analysis of flux predictions
fullReport.rxnID = [rxnTargetUptList; rxnTargetSecList];
fullReport.mean = [rxnTargetUpt_mean; rxnTargetSec_mean];
fullReport.SD = [rxnTargetUpt_SD; rxnTargetSec_SD];

target = sign(fullReport.mean);
%mean experimental uptake rates below analytical chemistry precision limit
%considered zero
target(abs(fullReport.mean) < param.boundPrecisionLimit) = 0;

fullReport.target = target;

fluxRxn = predictedFlux;
fullReport.v = fluxRxn(findRxnIDs(model, fullReport.rxnID));
predict = sign(fullReport.v);
%if the predicted flux magnitude is below fluxEpsilon,
%the predicted flux is considered zero
predict(abs(predict) < param.boundPrecisionLimit) = 0;
fullReport.predict = predict;

fullReport.agree = fullReport.target == fullReport.predict;
fullReport.dv = abs(fullReport.mean - fullReport.v);

%% Flux prediction accuracy

nRows = 1;
varTypes = {'double', 'double'};
varNames = {'wEuclidNorm', 'accuracy'};
comparisonStats = table('Size', [nRows length(varTypes)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%v = fullReport.v;
dv = fullReport.dv;
vExp = fullReport.mean;
bool = ~(isoutlier(dv) | isnan(dv) | isnan(vExp));
w  = sparse(1./(1 + (vExp.^2)));
if any(bool)
    comparisonStats.wEuclidNorm = sqrt(dv(bool)' * diag(w(bool)) * dv(bool));
else
    comparisonStats.wEuclidNorm = NaN;
end
% confusionMatrix
C = confusionmat(fullReport.target, fullReport.predict, 'ORDER',[1, 0, -1]);
comparisonStats.accuracy = sum(diag(C), 1) / sum(sum(C, 1));

if isnan(comparisonStats.wEuclidNorm)
    pause(0.1)
end
comparisonData.fullReport = fullReport;
comparisonData.comparisonStats = comparisonStats;
end