% running perturbations (complex I and V inhibition, media change to galactose
% from glucose) for iDopaNeuro1 and iDopaNeuro1.1 with all
% metabolomics constraints from control included; perturbed reaction's flux
% is gradually reduced and model is optimised at each step; metabolites
% that respond to the flux reduction can be analysed

% Set solver
clear
initCobraToolbox(0)
[~, ~] = changeCobraSolver('mosek', 'all', 0);

% load models
load('~/work/sbgCloud/exoMetDN/results/codeResults/iDN1/iDopaNeuro1Core/iDopaNeuro1Core.mat', 'iDopaNeuro1Core');
load('~/work/sbgCloud/exoMetDN/results/codeResults/iDN1/iDopaNeuro1ConditionType/iDopaNeuro1ConditionType.mat');

% using unWeightedTCBMfluxConc objective
iDN1 = iDopaNeuro1Core;
iDN1 = changeRxnBounds(iDN1, 'EX_gal[e]', 0, 'b'); %galactose not present

if isfield(iDN1, 'g0')
    iDN1 = rmfield(iDN1, 'g0');
end
if isfield(iDN1, 'g1')
    iDN1 = rmfield(iDN1, 'g1');
end

iDN11 = iDopaNeuro1ConditionType;
iDN11 = changeRxnBounds(iDN11, 'EX_gal[e]', 0, 'b'); %galactose not present

if isfield(iDN11, 'g0')
    iDN11 = rmfield(iDN11, 'g0');
end
if isfield(iDN11, 'g1')
    iDN11 = rmfield(iDN11, 'g1');
end

tcbmParam.method = 'fluxConc';
tcbmParam.solver = 'mosek';
tcbmParam.printLevel = 0;

%% complex V inhibition

model = iDN1;

model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;
        
iDN1_solutions.ATPS4m = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
iDN1_solutions.ATPS4m.Properties.VariableNames = {'Rxns', 'rxnFormulas', 'lb', 'ub'};

% original flux (max)
[solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
max_iDN1 = solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'ATPS4m'));
% minimum feasible flux (min) - currently not working
%model_min = changeObjective(model, 'ATPS4m');

%[solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
%min_iDN1 = solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'ATPS4m'));

for i = max_iDN1:-(max_iDN1/100):0
    try
        
        model = changeRxnBounds(model, 'ATPS4m', i, 'u');
        flux = i;
        [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
        iDN1_solutions.ATPS4m(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
        iDN1_solutions.ATPS4m.Properties.VariableNames{end} = ['ATPS4m flux = ' num2str(flux)];
        
    catch
        disp('--------------------------------------------------------------')
        disp(' ')
        disp(['iDN1 model infeasible at ATPS4m flux (TCBM) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')      
        sol_min = optimizeCbModel(changeObjective(iDN1, 'ATPS4m'), 'min');
        flux = sol_min.f;
        disp(['iDN1 model minimum feasible flux through ATPS4m (FBA) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')
    end
end

% iDN1.1 model
model = iDN11;

model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;
        
iDN11_solutions.ATPS4m = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
iDN11_solutions.ATPS4m.Properties.VariableNames = {'Rxns', 'rxnFormulas', 'lb', 'ub'};

% original flux (max)
[solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
max_iDN11 = solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'ATPS4m'));

for i = max_iDN11:-(max_iDN11/100):0
    try
        
        model = changeRxnBounds(model, 'ATPS4m', i, 'u');
        flux = i;
        [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
        iDN11_solutions.ATPS4m(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
        iDN11_solutions.ATPS4m.Properties.VariableNames{end} = ['ATPS4m flux = ' num2str(flux)];
        
    catch
        disp('--------------------------------------------------------------')
        disp(' ')
        disp(['iDN1.1 model infeasible at ATPS4m flux (TCBM) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')      
        sol_min = optimizeCbModel(changeObjective(iDN11, 'ATPS4m'), 'min');
        flux = sol_min.f;
        disp(['iDN1.1 model minimum feasible flux through ATPS4m (FBA) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')
    end
end

%% complex I inhibition (NADH2_u10m)

model = iDN1;

model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;
        
iDN1_solutions.NADH2_u10m = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
iDN1_solutions.NADH2_u10m.Properties.VariableNames = {'Rxns', 'rxnFormulas', 'lb', 'ub'};

% original flux (max)
[solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
max_iDN1 = solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'NADH2_u10m'));

for i = max_iDN1:-(max_iDN1/100):0
    try
        
        model = changeRxnBounds(model, 'NADH2_u10m', i, 'u');
        flux = i;
        [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
        iDN1_solutions.NADH2_u10m(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
        iDN1_solutions.NADH2_u10m.Properties.VariableNames{end} = ['NADH2_u10m flux = ' num2str(flux)];
        
    catch
        disp('--------------------------------------------------------------')
        disp(' ')
        disp(['iDN1 model infeasible at NADH2_u10m flux (TCBM) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')      
        sol_min = optimizeCbModel(changeObjective(iDN1, 'NADH2_u10m'), 'min');
        flux = sol_min.f;
        disp(['iDN1 model minimum feasible flux through NADH2_u10m (FBA) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')
    end
end

% iDN1.1 model
model = iDN11;

model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;
        
iDN11_solutions.NADH2_u10m = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
iDN11_solutions.NADH2_u10m.Properties.VariableNames = {'Rxns', 'rxnFormulas', 'lb', 'ub'};

% original flux (max)
[solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
max_iDN11 = solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'NADH2_u10m'));

for i = max_iDN11:-(max_iDN11/100):0
    try
        
        model = changeRxnBounds(model, 'NADH2_u10m', i, 'u');
        flux = i;
        [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
        iDN11_solutions.NADH2_u10m(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
        iDN11_solutions.NADH2_u10m.Properties.VariableNames{end} = ['NADH2_u10m flux = ' num2str(flux)];
        
    catch
        disp('--------------------------------------------------------------')
        disp(' ')
        disp(['iDN1.1 model infeasible at NADH2_u10m flux (TCBM) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')      
        sol_min = optimizeCbModel(changeObjective(iDN11, 'NADH2_u10m'), 'min');
        flux = sol_min.f;
        disp(['iDN1.1 model minimum feasible flux through NADH2_u10m (FBA) = ' num2str(flux) ' umol/gDW/h'])
        disp(' ')
    end
end

%% galactose media

% first iDN1 model
% media with glucose and no galactose

model = iDN1;
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present

model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;

iDN1_solutions.Gal = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
iDN1_solutions.Gal.Properties.VariableNames = {'Rxns', 'rxnFormulas', 'lb', 'ub'};

try
    [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
    iDN1_solutions.Gal(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
    iDN1_solutions.Gal.Properties.VariableNames{end} = 'Glucose media';
    disp('--------------------------------------------------------------')
    disp('Glucose media:')
    disp(['iDN1 model glucose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_glc_D[e]'))) ...
        ' umol/gDW/h'])
    disp(['iDN1 model galactose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_gal[e]'))) ...
        ' umol/gDW/h'])
    disp(' ')
catch
    disp('--------------------------------------------------------------')
    disp('Glucose media:')
    disp(['iDN1 model infeasible with glucose flux range between ' ...
        num2str(model.lb(findRxnIDs(model, 'EX_glc_D[e]'))) ' and ' ...
        num2str(model.ub(findRxnIDs(model, 'EX_glc_D[e]'))) ' umol/gDW/h'])
    disp(' ')
end

% glucose not present
model.lb(findRxnIDs(model, 'EX_gal[e]')) = model.lb(findRxnIDs(model, 'EX_glc_D[e]'));
model.ub(findRxnIDs(model, 'EX_gal[e]')) = model.ub(findRxnIDs(model, 'EX_glc_D[e]'));
model = changeRxnBounds(model, 'EX_glc_D[e]', 0, 'b');

try
    [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
    iDN1_solutions.Gal(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
    iDN1_solutions.Gal.Properties.VariableNames{end} = 'Galactose media';
    disp('--------------------------------------------------------------')
    disp('Galactose media:')
    disp(['iDN1 model glucose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_glc_D[e]'))) ...
        ' umol/gDW/h'])
    disp(['iDN1 model galactose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_gal[e]'))) ...
        ' umol/gDW/h'])
    disp(' ')
catch
    disp('--------------------------------------------------------------')
    disp('Galactose media:')
    disp(['iDN1 model infeasible with galactose flux range between ' ...
        num2str(model.lb(findRxnIDs(model, 'EX_gal[e]'))) ' and ' ...
        num2str(model.ub(findRxnIDs(model, 'EX_gal[e]'))) ' umol/gDW/h'])
    disp(' ')
end

% do the same for model iDN1.1
model = iDN11;
model = changeRxnBounds(model, 'EX_gal[e]', 0, 'b'); %galactose not present

model.osenseStr = 'min';
model.cf = 0;
model.cr = 0;
model.g = 2;
model.u0 = 0;
model.f = 1;

iDN11_solutions.Gal = table(model.rxns, model.rxnFormulas, model.lb, model.ub);
iDN11_solutions.Gal.Properties.VariableNames = {'Rxns', 'rxnFormulas', 'lb', 'ub'};

try
    [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
    iDN11_solutions.Gal(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
    iDN11_solutions.Gal.Properties.VariableNames{end} = 'Glucose media';
    disp('--------------------------------------------------------------')
    disp('Glucose media:')
    disp(['iDN1.1 model glucose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_glc_D[e]'))) ...
        ' umol/gDW/h'])
    disp(['iDN1.1 model galactose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_gal[e]'))) ...
        ' umol/gDW/h'])
    disp(' ')
catch
    disp('--------------------------------------------------------------')
    disp('Glucose media:')
    disp(['iDN1.1 model infeasible with glucose flux range between ' ...
        num2str(model.lb(findRxnIDs(model, 'EX_glc_D[e]'))) ' and ' ...
        num2str(model.ub(findRxnIDs(model, 'EX_glc_D[e]'))) ' umol/gDW/h'])
    disp(' ')
end

%glucose not present
model.lb(findRxnIDs(model, 'EX_gal[e]')) = model.lb(findRxnIDs(model, 'EX_glc_D[e]'));
model.ub(findRxnIDs(model, 'EX_gal[e]')) = model.ub(findRxnIDs(model, 'EX_glc_D[e]'));
model = changeRxnBounds(model, 'EX_glc_D[e]', 0, 'b');

try
    [solutions.unWeightedTCBMfluxConc, ~] = entropicFluxBalanceAnalysis(model, tcbmParam);
    iDN11_solutions.Gal(:,end+1) = table(solutions.unWeightedTCBMfluxConc.v);
    iDN11_solutions.Gal.Properties.VariableNames{end} = 'Galactose media';
    disp('--------------------------------------------------------------')
    disp('Galactose media:')
    disp(['iDN1.1 model glucose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_glc_D[e]'))) ...
        ' umol/gDW/h'])
    disp(['iDN1.1 model galactose flux = ' ...
        num2str(solutions.unWeightedTCBMfluxConc.v(findRxnIDs(model, 'EX_gal[e]'))) ...
        ' umol/gDW/h'])
    disp(' ')
catch
    disp('--------------------------------------------------------------')
    disp('Galactose media:')
    disp(['iDN1.1 model infeasible with galactose flux range between ' ...
        num2str(model.lb(findRxnIDs(model, 'EX_gal[e]'))) ' and ' ...
        num2str(model.ub(findRxnIDs(model, 'EX_gal[e]'))) ' umol/gDW/h'])
    disp(' ')
end

%% summarise results and calculate covariance

[bool,ExIDs] = ismember(iDN1.rxns(findExcRxns(iDN1)), iDN1.rxns);

IDs = ExIDs(bool);

for i=1:length(iDN1.rxns)
    [R,p] = corrcoef(table2array(iDN1_solutions.ATPS4m(i,5:end))', ...
    table2array(iDN1_solutions.ATPS4m(findRxnIDs(iDN1, 'ATPS4m'),5:end))');
    r_c5(i,1) = R(2,1);
    p_c5(i,1) = p(2,1);
end

max_dC5 = table2array(iDN1_solutions.ATPS4m(:,end))./table2array(iDN1_solutions.ATPS4m(:,5));
complex5_exRxns_summary_iDN1 = table(iDN1.rxns, iDN1.rxnNames, r_c5, table2array(iDN1_solutions.ATPS4m(:,5)), table2array(iDN1_solutions.ATPS4m(:,end)), max_dC5);
complex5_exRxns_summary_iDN1.Properties.VariableNames = {'RxnID', 'RxnNames', 'r (correlation with ATPS4m)', 'flux at max ATPS4m', 'flux at min ATPS4m', 'ratio min/max flux'};
complex5_sig_iDN1 = complex5_exRxns_summary_iDN1(find(p_c5<0.05),:);


for i=1:length(iDN1.rxns)
    [R,p] = corrcoef(table2array(iDN1_solutions.NADH2_u10m(i,5:end))', ...
    table2array(iDN1_solutions.NADH2_u10m(findRxnIDs(iDN1, 'NADH2_u10m'),5:end))');
    r_c1(i,1) = R(2,1);
    p_c1(i,1) = p(2,1);
end

max_dC1 = table2array(iDN1_solutions.NADH2_u10m(:,end))./table2array(iDN1_solutions.NADH2_u10m(:,5));
complex1_exRxns_summary_iDN1 = table(iDN1.rxns, iDN1.rxnNames, r_c1, table2array(iDN1_solutions.NADH2_u10m(:,5)), table2array(iDN1_solutions.NADH2_u10m(:,end)), max_dC1, max_dC1.*sign(complex5_exRxns_summary_iDN1.("flux at min ATPS4m")));
complex1_exRxns_summary_iDN1.Properties.VariableNames = {'RxnID', 'RxnNames', 'r (correlation with NADH2_u10m)', 'flux at max NADH2_u10m', 'flux at min NADH2_u10m', 'ratio min/max flux', 'ratio min/max flux, flux sign'};
complex1_sig_iDN1 = complex5_exRxns_summary_iDN1(find(p_c1<0.05),:);

% low flux values excluded:

param.treshold_lb = -0.001;
param.treshold_ub = 0.001;
label_value = {};
for i=1:length(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m"))
    if sign(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i)) == ...
            sign(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i))
        if abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i)) > ...
                abs(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i))
            if round(abs(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i)) / ...
                    abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i)),3) == 1
                label_value{i,1} = ...
                    [num2str(abs(round(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i),3))) ...
                    ' (0%)'];
            else
                decrease = (abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i)) ...
                    - abs(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i))) / ...
                    abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i));
                label_value{i,1} = ...
                    [num2str(abs(round(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i),3))) ...
                    ' (-' num2str(round(decrease*100)) '%)'];
            end
        elseif abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i)) < ...
                abs(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i))
            if round(abs(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i)) / ...
                    abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i)),3) == 1
                label_value{i,1} = ...
                    [num2str(abs(round(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i),3))) ...
                    ' (0%)'];
            else
                increase = (abs(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i)) ...
                    - abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i))) / ...
                    abs(complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m")(i));
                label_value{i,1} = ...
                    [num2str(abs(round(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i),3))) ...
                    ' (+' num2str(round(increase*100)) '%)'];
            end
        else
            label_value{i,1} = ...
                [num2str(abs(round(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i),3))) ...
                ' (0%)'];
        end
    else
        label_value{i,1} = ...
            [num2str(abs(round(complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m")(i),3))) ...
            ' (sign change!)'];
    end
end
param.EdgeLabel.text = label_value;
param.EdgeLabel.rxnID = complex1_exRxns_summary_iDN1.RxnID;

[graph, summary] = metUtilisation(iDN1, 'fad[m]', complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m"),1, param);

atp_flux = (summary.flux_v.*summary.scoff);
atp_consumption = sum(atp_flux(find(atp_flux<0)));

param.treshold_lb = 0.05 * atp_consumption;
param.treshold_ub = -0.05 * atp_consumption;

[graph2, summary2] = metUtilisation(iDN1, 'fad[m]', complex1_exRxns_summary_iDN1.("flux at max NADH2_u10m"),1, param);

figure();
H = plot(graph2.graph,'NodeLabel',graph2.nLabels,'NodeFontSize',10,'EdgeColor',graph2.eColour, ...
    'EdgeLabel', graph2.edgeLabels,'LineWidth',graph2.LWidths, ...
    'ArrowSize', 10, 'ArrowPosition', 0.65, 'Interpreter','none');
box off
layout(H,'layered','Direction','right')

RxnNames1 = iDN1.rxnNames(findRxnIDs(iDN1, H.NodeLabel(H.XData == 1)));
RxnNames3 = iDN1.rxnNames(findRxnIDs(iDN1, H.NodeLabel(H.XData == 3)));
text(H.XData(H.XData == 2), H.YData(H.XData == 2)+0.5, H.NodeLabel(H.XData == 2), ...
    'VerticalAlignment','bottom', 'HorizontalAlignment', 'center', 'FontSize', 15, ...
    'Interpreter','none')
text(H.XData(H.XData == 1)-0.05, H.YData(H.XData == 1), RxnNames1, ...
    'VerticalAlignment','middle', 'HorizontalAlignment', 'right', 'FontSize', 10, ...
    'Interpreter','none')
text(H.XData(H.XData == 3)+0.05, H.YData(H.XData == 3), RxnNames3, ...
    'VerticalAlignment','middle', 'HorizontalAlignment', 'left', 'FontSize', 10, ...
    'Interpreter','none')
H.NodeLabel = {};
line(NaN,NaN,'Color','#A2142F','LineStyle','-')
line(NaN,NaN,'Color','#0072BD','LineStyle','-')
%legend([metID ' utilisation network:'], 'production', 'consumption')

%% plot specific intracellular changes (barplot)

Rxns = {'NADH2_u10m', 'FADH2ETC', 'ATPS4m','PYK', 'PGL'};
Labels = {'Complex I', 'Complex II', 'Complex V', 'Pyruvate kinase', '6-Phosphogluconolactonase'};
fluxData = complex1_exRxns_summary_iDN1;
[LIA,LOCB] = ismember(Rxns,fluxData.RxnID);
selectedData = fluxData(LOCB,:);

figure
bar(1:5, [abs(selectedData.("flux at max NADH2_u10m")) abs(selectedData.("flux at min NADH2_u10m"))], 1)

ylim([0 55])
box off
set(gca, 'XTick', 1:5, 'XTickLabel', Labels)
% Add title and axis labels
%title('Predicted flux of the core energy-related reactions in control and complex I inhibition.')
ylabel('Flux (uMol/gDW/hr)')
% Add a legend
legend('Control', 'Complex I inhibition')
%% to be continued: calculate accuracies and plot sigmoidal figure

exoMetabolomicsDataDir = ['~' filesep 'work' filesep 'sbgCloud' filesep ...
    'exoMetDN' filesep 'data' filesep 'omics'];

% Read the files
exometData = readtable([exoMetabolomicsDataDir filesep ...
    '20210118_validationData_K7_final_ordered_umol_gDW_h.csv']);
perturbationData = readtable([exoMetabolomicsDataDir filesep ...
    '20210118_validationData_K7_finalAll_ordered_umol_gDW_h.csv']);
rxns = cell(size(perturbationData.variable));
[mlt, nlt] = size(exometData.variable);

% Complete table
for i = 1:mlt
    rxns(ismember(perturbationData.variable, exometData.variable{i})) = ...
        exometData.exRxns(i);
end
perturbationData = [perturbationData rxns];
perturbationData.Properties.VariableNames{'Var10'} = 'rxnID';
perturbationData.Properties.VariableNames{'sds'} = 'SD';
perturbationData.Properties.VariableNames{'variable'} = 'name';
perturbationData.Properties.VariableNames{'platform'} = 'Platform';

% c1_validationData
rotenoneBool = strcmp(perturbationData.condition, 'rotenone');
C1ValidationData = perturbationData(rotenoneBool, :);
C1ValidationData = sortrows(C1ValidationData, 'rxnID');
% c5_validationData
rotenoneBool = strcmp(perturbationData.condition, 'oligomycin');
C5ValidationData = perturbationData(rotenoneBool, :);
C5ValidationData = sortrows(C5ValidationData, 'rxnID');
% galactose_validationData
galactoseBool = strcmp(perturbationData.condition, 'galactose');
galValidationData = perturbationData(galactoseBool, :);
galValidationData = sortrows(galValidationData, 'rxnID');

% glucose_validationData
glucoseBool = strcmp(perturbationData.condition, 'control');
glcValidationData = perturbationData(glucoseBool, :);
glcValidationData = sortrows(glcValidationData, 'rxnID');

predictedFlux = complex1_exRxns_summary_iDN1.("flux at min NADH2_u10m");
validationData = C1ValidationData;
comparisonData_C1 = quanQualAcc(iDN1, predictedFlux, validationData);

predictedFlux = complex5_exRxns_summary_iDN1.("flux at min ATPS4m");
validationData = C5ValidationData;
comparisonData_C5 = quanQualAcc(iDN1, predictedFlux, validationData);

predictedFlux = iDN1_solutions.Gal.("Galactose media");
validationData = galValidationData;
comparisonData_gal = quanQualAcc(iDN1, predictedFlux, validationData);

predictedFlux = iDN1_solutions.Gal.("Glucose media");
validationData = glcValidationData;
comparisonData_glc = quanQualAcc(iDN1, predictedFlux, validationData);

summaryTableAccuracies = [{"Control"}, comparisonData_glc.comparisonStats; ...
    {"Galactose"}, comparisonData_gal.comparisonStats;...
    {"Complex I inhibition"}, comparisonData_C1.comparisonStats;...
    {"Complex V inhibition"}, comparisonData_C5.comparisonStats];

% plot predicted and measured exchanges
exoMet = glcValidationData;
fullReport = comparisonData_C1.fullReport;
condition = 'Complex I inhibition';
objective = 'unWeightedTCBMfluxConc';
comparisonObjective=[];
labelType = 'metabolite'; %'platform'
saveFigures = 1;
comparison = []; % second predicted flux to be added for comparison
comparison_label = 'control';
driver_plotSigmoidalFigure

%calculate and plot delta of the changes
measured_control_data = glcValidationData;
measured_test_data = C1ValidationData;
predicted_control_data = comparisonData_glc.fullReport;
predicted_test_data = comparisonData_C1.fullReport;
driver_perturbationsDeltaChange

exoMet = C1ValidationData;
exoMet.mean = measured_dv;
exoMet.SD = measured_dvsd;
exoMet.sign = measured_label;
fullReport = comparisonData_C1.fullReport;
fullReport.mean = predicted_dv;
fullReport.sign = predicted_label;
condition = 'complex I inhibition';
objective = 'unWeightedTCBMfluxConc';
comparisonObjective=[];
labelType = 'metabolite'; %'platform'
saveFigures = 1;

driver_plotSigmoidalFigure_deltaChange


