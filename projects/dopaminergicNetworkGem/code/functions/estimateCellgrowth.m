function [cellConcentration, timing] = estimateCellgrowth(extraCellParams, objFunCell) 
%% Cell growth estimation 
% estimate the growth of a cell pupulation based of measurements at
% different time points 

% Averina Nicolae 2015. LCSB, Belval, Luxembourg. 
% Jennifer Modamio code addaptation 03.2018, LCSB, Belval,Luxembourg. 

% INPUTS 
%
%    structure type: 
%    Example: extraCellParams
% 
%    extraCellParams.cellModel = string type growth 
%                                'linear' // 'constant' // 'exponential' 
%
%    extraCellParams.tsampled  = cell containing a double with timing;
%                                [0;12;18] (days)
%
%    extraCellParams.cellsampled = double with cell density in each
%    timepoint 
%                                [400000;763298;855355]
%
%    extraCellParams.cellSD = duble with number of seeded cells 
%                                [400000]
% 

% extraCellParams.cellModel = 'linear';
% extraCellParams.tsampled  = [0;12;18] ;
% extraCellParams.cellsampled =  [400000;763298;855355];
% extraCellParams.cellSD = [80000];
% extraCellParams.tapprox = [23];


% corrSDOption = 'min'; % correct SD values found to = 0
% changeCobraSolver('matlab', 'NLP')
% 
    

%% Parameter 
              
fminconopt = optimset('MaxFunEvals', 1e7, 'TolFun', 1e-8, 'display', 'off'); 
objFunCell = 'cellFitSQD'; 

cellp0 = [extraCellParams.cellsampled(1), 0];
cellpmin = [extraCellParams.cellsampled(1) - extraCellParams.cellSD(1), 0];
cellpmax = [extraCellParams.cellsampled(1) + extraCellParams.cellSD(1), 0.5*max(extraCellParams.cellsampled)];
growthmodel = extraCellParams.cellModel ;

% extraCellParams.tsampled  = tsampleCells;
% extraCellParams.cellsampled = cellDensTotal;
% extraCellParams.cellSD = cellSDTotal;
% extraCellParams.tapprox = tapprox

%% solving problem 
[estimCellParam, fvalCell] = solveMetabolomicsFitting(objFunCell, cellp0, cellpmin, cellpmax, extraCellParams, fminconopt);

%% adapting to cell specific data 

tcellsim = [extraCellParams.tsampled(1):0.1:extraCellParams.tapprox];
[cellsim, tcellsim] = modelCellGrowth(growthmodel, estimCellParam, tcellsim);

cellConcentration = cellsim 
timing = tcellsim
end

%% 