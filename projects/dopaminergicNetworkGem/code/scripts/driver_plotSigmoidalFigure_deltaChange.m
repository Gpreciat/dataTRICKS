%% Plot change in uptake and secretion rates as well as error bars on logarithmic scale
%options to be specified and loaded before:
% exoMet = measured_delta;
% fullReport = predicted_delta;
% condition = 'Complex I inhibition';
% objective = 'unWeightedTCBMfluxConc';
% comparisonObjective=[];
% labelType = 'metabolite'; %'platform'
% saveFigures = 1;

close all

labelsL = strtrim(exoMet.name);
labelsR = exoMet.Platform;
       
nlt=length(labelsL);
bool=true(nlt,1);
maxCharacters=0;
for i=1:nlt
    nCharacters=length(labelsL{i});
    if nCharacters>maxCharacters
        maxCharacters=nCharacters;
    end
    if isempty(labelsL{i})
        labelsL{i}='';
        bool(i)=0;
    else
        if ismember(labelsL{i},{'Nc20:4 ';'Ex_icit[e]';'Nc20:4 ';...
                'EX_CE1557[e]';'EX_HC00900[e]';...
                'Transport of (R)-3-Hydroxybutanoate  via H+ Symport';...
                'EX_M03117[e]'; 'Proline'; 'Threonine'; 'Leucine'; ...
                'Lysine'; 'Methionine'; 'Isoleucine'})
            bool(i)=0;
        end
    end
    if isempty(labelsR{i})
        labelsR{i}='';
    end
end
labelsL = strrep(labelsL,'Exchange of ','');
labelsL = strrep(labelsL,'3, 4-dihydroxy','3,4-dihydroxy');
labelsL = strrep(labelsL,'N-Acetyl-Tyrosine','N-Acetyl-L-tyrosine');
labelsL = strrep(labelsL,'.','-');
labelsL = strrep(labelsL, 'Gamma-aminobutyric-acid', 'GABA');
labelsL = strrep(labelsL,'-acid',' acid');
labelsL = strrep(labelsL,'-Acid',' acid');
labelsL = strtrim(labelsL);
exoMet.labelsL=labelsL;

labelsR = strrep(labelsR,'GC-MS','$GCMS$');
labelsR = strrep(labelsR,'LC-MS','$AccQtag$');
labelsR = strrep(labelsR,'dmpa_old','$DmPa^{13}C$');
labelsR = strrep(labelsR,'dmpa_new','$DmPa^{2}H$');
labelsR = strrep(labelsR,'dmpaN','$DmPa^{2}H$');
labelsR = strrep(labelsR,'dmpa','$DmPa^{13}C$');
labelsR = strrep(labelsR,'bzcl_old','$BzCl$');
labelsR = strrep(labelsR,'bzcl_new','$BzCl$');
labelsR = strrep(labelsR,'bzcl','$BzCl$');
exoMet.labelsR=labelsR;


switch labelType
    case 'metabolitePlatform'
        for i=1:nlt
            numberOfCharacters = maxCharacters - length(labelsL{i}) + 1;
            numberOfCharacters = 1;
            labels{i}=[labelsL{i} ' ' repmat('\ ',1,numberOfCharacters) ' ' labelsR{i}];
        end
    case 'metabolite'
        labels = labelsL;
    case 'platform'
        labels = labelsR;
end
exoMet.labels=labels;

% Exchange rates, excluding certain metabolites
exoMet = exoMet(bool,:);
%sort by mean measured exchange rate
exoMet.meanOld = exoMet.mean;
exoMet.mean(contains(exoMet.sign, 'UU')) = -exoMet.mean(contains(exoMet.sign, 'UU'));
exoMet.mean(contains(exoMet.sign, '0U')) = -exoMet.mean(contains(exoMet.sign, '0U'));
exoMet.mean(contains(exoMet.sign, 'SU')) = -exoMet.mean(contains(exoMet.sign, 'SU'));
exoMet.mean(contains(exoMet.sign, 'SS')) = exoMet.mean(contains(exoMet.sign, 'SS'));
exoMet.mean(contains(exoMet.sign, '0S')) = exoMet.mean(contains(exoMet.sign, '0S'));
exoMet.mean(contains(exoMet.sign, 'US')) = exoMet.mean(contains(exoMet.sign, 'US'));
exoMet_dirChange_bool = contains(exoMet.sign, 'US')|contains(exoMet.sign, 'SU');
exoMet = sortrows(exoMet,'mean','ascend');

fullReport.meanOld = fullReport.mean;
fullReport.mean(contains(fullReport.sign, 'UU')) = -fullReport.mean(contains(fullReport.sign, 'UU'));
fullReport.mean(contains(fullReport.sign, '0U')) = -fullReport.mean(contains(fullReport.sign, '0U'));
fullReport.mean(contains(fullReport.sign, 'SU')) = -fullReport.mean(contains(fullReport.sign, 'SU'));
fullReport.mean(contains(fullReport.sign, 'SS')) = fullReport.mean(contains(fullReport.sign, 'SS'));
fullReport.mean(contains(fullReport.sign, '0S')) = fullReport.mean(contains(fullReport.sign, '0S'));
fullReport.mean(contains(fullReport.sign, 'US')) = fullReport.mean(contains(fullReport.sign, 'US'));
prediction_dirChange_bool = contains(fullReport.sign, 'US')|contains(fullReport.sign, 'SU');

%plotType ='linearExchanges';
plotType = 'logExchanges';
switch plotType
    case 'linearExchanges'
       
        figure
        hold on
        ylim([0,length(exoMet.mean)+0.5])
        errorbar(exoMet.mean,1:length(exoMet.mean),2*exoMet.SD,'horizontal','*','LineWidth',4)
        %find the predicted upakes in the full report
        fullReportBool = strcmp(fullReport.data,data);
        predictedExchange.rxnID = fullReport.rxnID(fullReportBool);
        predictedExchange.v = fullReport.v(fullReportBool);
        %handle potentially missing reactions in the model
        nPredictedExchanges = length(predictedExchange.v);
        predictedExchange.rxnID{nPredictedExchanges+1} = 'Missmatch';
        predictedExchange.v(nPredictedExchanges+1) = NaN;
        [~,predictedExLOCB] = ismember(exoMet.rxnID,predictedExchange.rxnID);
        predictedExLOCB(predictedExLOCB==0)= nPredictedExchanges+1;
        plot(predictedExchange.v(predictedExLOCB),1:nnz(predictedExLOCB),'*r')
        hold off
        legend('measured','predicted','Location','northwest');
        title('Exchange rates')
        xlabel('uMol/gDW/hr');
        yticks(1:length(exoMet.mean))
        %set(gca,'XScale','log');
        yticklabels(exoMet.labels)
        ax = gca;
        ax.TickLabelInterpreter='none';
        set(gca,'FontSize',14)
        ax = gca;
        ax.TickLabelInterpreter='none';
        set(gca,'FontSize',14)
        grid on
        ax = gca;
        ax.GridColor = [0 .5 .5];
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.5;
        ax.Layer = 'top';
        
        % Exchange rates, excluding those with large measured values
        bool = abs(exoMet.mean) <10;
        exoMet = exoMet(bool,:);
        
        figure
        hold on
        ylim([0,length(exoMet.mean)+0.5])
        errorbar(exoMet.mean,1:length(exoMet.mean),2*exoMet.SD,'horizontal','*','LineWidth',4)
        %find the predicted upakes in the full report
        fullReportBool = strcmp(fullReport.data,data);
        predictedExchange.rxnID = fullReport.rxnID(fullReportBool);
        predictedExchange.v = fullReport.v(fullReportBool);
        %handle potentially missing reactions in the model
        nPredictedExchanges = length(predictedExchange.v);
        predictedExchange.rxnID{nPredictedExchanges+1} = 'Missmatch';
        predictedExchange.v(nPredictedExchanges+1) = NaN;
        [~,predictedExLOCB] = ismember(exoMet.rxnID,predictedExchange.rxnID);
        predictedExLOCB(predictedExLOCB==0)= nPredictedExchanges+1;
        plot(predictedExchange.v(predictedExLOCB),1:nnz(predictedExLOCB),'*r')
        hold off
        legend('measured','predicted','Location','northwest');
        title('Exchange rates')
        xlabel('uMol/gDW/hr');
        yticks(1:length(exoMet.mean))
        %set(gca,'XScale','log');
        yticklabels(exoMet.labels)
        ax = gca;
        ax.TickLabelInterpreter='none';
        set(gca,'FontSize',10)
        grid on
        ax = gca;
        ax.GridColor = [0 .5 .5];
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.5;
        ax.Layer = 'top';
        
    case 'logExchanges'
             
        %sort by the mean measured exchange for all except possibly zero
        %change
        zerobool = (exoMet.mean-exoMet.SD)<=0 & (exoMet.mean+exoMet.SD)>=0;
        posbool = sign(exoMet.mean)==1 & ~zerobool;
        negbool = sign(exoMet.mean)==-1 & ~zerobool;
        exoMet.uptakeProbability = normcdf(0,exoMet.mean,exoMet.SD);
        
        toRankOrder = zeros(length(exoMet.mean),1);
        toRankOrder(zerobool) = exoMet.uptakeProbability(zerobool);
        toRankOrder(posbool) = exoMet.mean(posbool)+1;
        toRankOrder(negbool)= exoMet.mean(negbool);
        
        [sortedUptakeProbability,xi] = sort(toRankOrder);
        exoMet = exoMet(xi,:);
        
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on       
        
        %errorbar(exoMet.mean,1:length(exoMet.mean),2*exoMet.SD,'horizontal','*','LineWidth',4)
        x = exoMet.mean;
        xl = exoMet.mean - exoMet.SD;
        xu = exoMet.mean + exoMet.SD;
        
        %log modulus transformation
        x = logmod(x,10);
        xl = logmod(xl,10);
        xl = x - xl;
        xu = logmod(xu,10);
        xu = xu - x;
       
        exoMet.x=x;
        exoMet.xl=xl;
        exoMet.xu=xu;
        
        x_dC = ones(length(exoMet.x),1)*NaN;
        x_dC(exoMet_dirChange_bool) = exoMet.x(exoMet_dirChange_bool);
        xl_dC = ones(length(exoMet.x),1)*NaN;
        xl_dC(exoMet_dirChange_bool) = exoMet.xl(exoMet_dirChange_bool);
        xu_dC = ones(length(exoMet.x),1)*NaN;
        xu_dC(exoMet_dirChange_bool) = exoMet.xu(exoMet_dirChange_bool);
        
        %find the predicted exchanges in the full report
        %fullReportBool = strcmp(fullReport.data,data) & strcmp(fullReport.objective,objective) & (strcmp(fullReport.model,'modelSec') | strcmp(fullReport.model,'modelUpt'));
        predictedExchange=table;
        predictedExchange.rxnID = fullReport.rxnID;
        predictedExchange.v = fullReport.v;
        predictedExchange.vlogmod = logmod(predictedExchange.v,10);
        
        %place the predictions in the same order as the measurements
        [LIA,LOCB] = ismember(predictedExchange.rxnID,exoMet.rxnID);
        %sanity check
        exoMet.rxnID2(LOCB(LIA)) =  predictedExchange.rxnID(LIA);
        v=ones(length(exoMet.x),1)*NaN;
        v(LOCB(LIA)) = logmod(predictedExchange.v(LIA),10);
        v_dC = ones(length(exoMet.x),1)*NaN;
        v_dC(LOCB(and(LIA, prediction_dirChange_bool))) = logmod(predictedExchange.v(and(LIA, prediction_dirChange_bool)),10);
        %log modulus transformation
        %v = logmod(v,10);
        exoMet.vPredicted = v;
        
        if ~isempty(comparison)
            %find the predicted exchanges in the full report
            predictedExchange2.rxnID = comparison.rxnID;
            predictedExchange2.v = fullReport.v;
            
            %place the predictions in the same order as the measurements
            [LIA,LOCB] = ismember(predictedExchange2.rxnID,exoMet.rxnID);
            %sanity check
            %exoMet.rxnID2(LOCB(LIA)) =  predictedExchange.rxnID(LIA);
            v2=ones(length(exoMet.x),1)*NaN;
            v2(LOCB(LIA)) = logmod(predictedExchange2.v(LIA),10);
            
            %log modulus transformation
            %v2 = logmod(v2,10);
            
            %limits
            min_v = min([v;v2;x-xl]);
            max_v = max([v;v2;x+xu]);
            bool=~isnan(v) & ~isnan(v2);
        else
            min_v = min([v;x-xl]);
            max_v = max([v;x+xu]);
            bool=~isnan(v);
        end

        %patch in the background
        X = [min_v max(x+xu) max(x+xu) min_v];

        Y = [nnz(negbool & bool)+nnz(zerobool & bool)+0.5 nnz(negbool & bool)+nnz(zerobool & bool)+0.5 length(exoMet.mean(bool))+0.5 length(exoMet.mean(bool))+0.5];
        patch(X,Y,[221/255, 192/255, 237/255],'FaceAlpha',0.5)
        
        Y = [nnz(negbool & bool)-0.5 nnz(negbool & bool)-0.5 nnz(negbool & bool)+nnz(zerobool & bool)+0.5 nnz(negbool & bool)+nnz(zerobool & bool)+0.5];
        patch(X,Y,[192/255, 237/255, 207/255],'FaceAlpha',0.5)
        
        Y = [0 0 nnz(negbool & bool)-0.5 nnz(negbool & bool)-0.5];
        patch(X,Y,[192/255, 233/255, 237/255],'FaceAlpha',0.9)
        
        
        
        %errorbar
        errorbar(x(bool),1:length(exoMet.mean(bool)),xl(bool),xu(bool),'horizontal','*','LineWidth',2)
        errorbar(x_dC(bool),1:length(exoMet.mean(bool)),xl_dC(bool),xu_dC(bool),'horizontal','*','LineWidth',2, 'Color', 'black')

        hold on
        
        %predicted fluxes
        plot(v(bool),1:length(v(bool)),'d','MarkerSize',12, 'MarkerEdgeColor', [0 0.4470 0.7410])%,'MarkerFaceColor',[1 1 1])
        %plot(v_dC(bool),1:length(v(bool)),'dk','MarkerSize',12)%,'MarkerFaceColor',[1 1 1])

        spearman = corr(x(isfinite(v)),v(isfinite(v)),'Type','Spearman');
        fprintf('%s%g\n',['Objective ' objective ', Spearman correlation = '],spearman);
%         spearman = corr(x(isfinite(v) & negbool),v(isfinite(v) & negbool),'Type','Spearman');
%         fprintf('%s%g\n',['Objective ' objective ' uptakes, Spearman correlation = '],spearman);
%         spearman = corr(x(isfinite(v) & posbool),v(isfinite(v) & posbool),'Type','Spearman');
%         fprintf('%s%g\n',['Objective ' objective ' secretions, Spearman correlation = '],spearman);
                
        %x & y limits

        ylim([0,length(exoMet.mean(bool))+0.5])
                xlabel('Change in flux $\delta$v (uMol/gDW/hr) on $\log_{10}$ scale: $\textrm{sign}(x) \times \log_{10}(1+x)$','interpreter','latex');

        %second predicted fluxes
        if ~isempty(comparison)                   
            %predicted fluxes
            plot(v2(bool),1:length(v2(bool)),'ok','MarkerSize',10)
            title(['Comparison of measured and predicted net exchange reaction rates. \color[rgb]{0 0.4470 0.7410}(' condition ') \color{black}(' comparison_label ')'])
            xlim([min_v max_v])
            %xlim([-0.5 0.5])
            legend({'measured secretion','no change', 'measured uptake','\delta v (measurement)','\delta v (prediction)',comparison_label,},'location','northwest');

        else
            xlim([min_v max_v])
            %'\delta v (prediction, sign change)'
            legend({'measured secretion (perturbation)','no change','measured uptake (perturbation)','\delta v (measurement)','\delta v (measurement, sign change)','\delta v (prediction)',},'location','northwest');
            %xlim([-0.5 0.5])
            title(['Comparison of measured and predicted changes in reaction rates, \delta v, between control and ' condition ')'])
        end
        
        %labels
        ax = gca;
        
        %xtick=[-100:10:100];
        %xticklab = cellstr(num2str(sign(xtick).*log10(1+abs(xtick)), '%d'));
        %set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex')

        ax.YTick = 1:length(exoMet.mean(bool));
        ax.YTickLabel = exoMet.labels(bool);
        ax.TickLabelInterpreter='latex';
        ax.FontSize = 14;
        switch labelType
            case 'metabolitePlatform'
                ax.YAxisLocation = 'right';
            case 'metabolite'
                ax.YAxisLocation = 'right';
            case 'platform'
                ax.YAxisLocation = 'left';
        end
        
%         %right labels
%         plot(v,1:nnz(predictedExLOCB),'*r')
%         axL = gca;
%         axL.YAxisLocation = 'right';
%         axL.Color = 'none';
%         axL.Box = 'off';
%         yticks(1:length(exoMet.mean))
%         axL.TickLabelInterpreter='latex';
%         axL.YColor ='black';
%         axL.FontSize = 10;
%         yt = yticks;
%         yticklabels(exoMet.labelsR)
        
        grid on
        ax = gca;
        ax.GridColor = [0 .5 .5];
        ax.GridLineStyle = '--';
        ax.GridAlpha = 0.5;
        ax.Layer = 'top';
 
        if 0
            xlim([-0.25 0.15])
            ylim([15.5 69.5])
            zoomed = 'zoom';
        else
            zoomed = '';
        end

        if saveFigures
            savefig(['iDN1_measured_vs_predicted_exchanges_' regexprep(condition, ' +', '_') '_' labelType zoomed])
            saveas(gcf,['iDN1_measured_vs_predicted_exchanges_' regexprep(condition, ' +', '_') '_' labelType zoomed],'png')
        end

    case 'logUptSec'
        log10AbsMean = log10(abs(exoMet.mean));
        posbool = sign(exoMet.mean)==1;
        negbool = sign(exoMet.mean)==-1;
        log10AbsSD = log10(exoMet.SD);
        
        % Uptake rates, excluding certain metabolites
        negbool = negbool & bool;
        figure
        hold on
        errorbar(-exoMet.mean(negbool),1:nnz(negbool),exoMet.SD(negbool),'horizontal','*','LineWidth',4)
        %find the predicted upakes in the full report
        fullReportBool = strcmp(fullReport.data,data) & strcmp(fullReport.objective,objective) & strcmp(fullReport.model,'modelSec');
        predictedUptake.rxnID = fullReport.rxnID(fullReportBool);
        predictedUptake.v = fullReport.v(fullReportBool);
        %handle potentially missing reactions in the model
        nPredictedUptakes = length(predictedUptake.v);
        predictedUptake.rxnID{nPredictedUptakes+1} = 'Missmatch';
        predictedUptake.v(nPredictedUptakes+1) = NaN;
        [~,predictedNegLOCB] = ismember(exoMet.rxnID(negbool),predictedUptake.rxnID);
        predictedNegLOCB(predictedNegLOCB==0)= nPredictedUptakes+1;
        plot(-predictedUptake.v(predictedNegLOCB),1:nnz(negbool),'*r')
        hold off
        legend('measured','predicted','Location','southwest');
        title('Uptake rates')
        xlabel('log(uMol/gDW/hr)');
        yticks(1:nnz(negbool))
        set(gca,'XScale','log');
        yticklabels(exoMet.labels(negbool))
        ax = gca;
        ax.TickLabelInterpreter='none';
        set(gca,'FontSize',10)
        %%
        
        % Secretion rates, excluding certain metabolites
        posbool = posbool & bool;
        figure
        hold on
        errorbar(exoMet.mean(posbool),1:nnz(posbool),exoMet.SD(posbool),'horizontal','*','LineWidth',4)
        ylim([0,nnz(posbool)+0.5])
        %find the predicted secretions in the full report
        fullReportBool = strcmp(fullReport.data,data) & strcmp(fullReport.objective,objective) & strcmp(fullReport.model,'modelUpt');
        predictedSecretion.rxnID = fullReport.rxnID(fullReportBool);
        predictedSecretion.v = fullReport.v(fullReportBool);
        %handle potentially missing reactions in the model
        nPredictedSecretions = length(predictedSecretion.v);
        predictedSecretion.rxnID{nPredictedSecretions+1} = 'Missmatch';
        predictedSecretion.v(nPredictedSecretions+1) = NaN;
        [~,predictedPosLOCB] = ismember(exoMet.rxnID(posbool),predictedSecretion.rxnID);
        predictedPosLOCB(predictedPosLOCB==0)= nPredictedSecretions+1;
        plot(-predictedSecretion.v(predictedPosLOCB),1:nnz(posbool),'*r')
        hold off
        legend('measured','predicted','Location','northeast');
        title('Secretion rates')
        xlabel('log(uMol/gDW/hr)');
        yticks(1:nnz(posbool))
        yticklabels(exoMet.labels(posbool))
        ax = gca;
        ax.TickLabelInterpreter='none';
        set(gca,'XScale','log');
        set(gca,'FontSize',10)
end