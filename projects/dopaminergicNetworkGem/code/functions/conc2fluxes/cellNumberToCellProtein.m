%this is the raw data obtained from Can
dataFolder = '~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/2019_conc2fluxes/';

fileName = '201122_metabolomicsData_Cornelius_vs_Can.xlsx';

close all

[~, ~, X] = xlsread([dataFolder fileName]);

[mlt,nlt]=size(X);

X{1,5}='mean';
for i=2:mlt
    X{i,5} = (X{i,2} +X{i,3})/2;
    corn(i)=X{i,5};
end
corn=corn(2:end);

X{1,10}='mean';
for i=2:mlt
    X{i,10} = (X{i,7} +X{i,8})/2;
    can(i)=X{i,10};
end
can=can(2:end);

X{1,11}='can2corn';
for i=2:mlt
    X{i,11} = X{i,5}/X{i,10};
    can2corn(i)= X{i,5}/X{i,10};
end
can2corn=can2corn(2:end);
can2CornOutlier = isoutlier(can2corn);

figure
plot(can2corn,[1:length(can2corn)],'.')
yticks([1:length(can2corn)])
set(gca,'TickLabelInterpreter','none')
labels=X(2:end,1);
yticklabels(labels)

can2corn(can2CornOutlier);

factorScale = mean(can2corn(~can2CornOutlier))
% factorScale = 0.142701415473448
%this approximately scales cans data to the scale of cornelius' data
X{1,12}='canCorned';
for i=2:mlt
    X{i,12} = X{i,10}*factorScale;
end

figure
hold on
ylt=length(corn(~can2CornOutlier));
plot(corn(~can2CornOutlier),[1:ylt],'*')
plot(can(~can2CornOutlier),[1:ylt],'*r')
plot(can(~can2CornOutlier)*factorScale,[1:ylt],'*g')
yticks([1:ylt])
ylim([0,ylt+1])
xlim([-20,20])
legend({'protein abundance normalisation','cell count normalisation','cell count data, renormalised to protein abundance'})
set(gca,'TickLabelInterpreter','none')
yticklabels(labels(~can2CornOutlier))
saveas(gcf,['~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/2019_conc2fluxes/' datestr(now,30) '_renormalised_metabolomics'],'jpg')