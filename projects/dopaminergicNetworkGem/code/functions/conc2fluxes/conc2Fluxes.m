% Code to convert extracellular concentration changes to exchange fluxes
% Differentiation protocol from NESC (neuroepithelial stem cells) to DN
% (dopaminergic neurons).
%
% The cell count and metabolic phenotype of the cell population changes over time.
% We only consider changes in concentration from day 19 to 23 as the cell population
% number reaches a plateau.

%Adapted from
%~/work/sbgCloud/programReconstruction/projects/exoMetDN/code/exometabolomics/conc2fluxes_mathmodel_2018/exometabolomicFluxStimation.m

% List of samples provided by Edinson
% 'Exo-metabolomics sample list_26-1-15.xlsx'

%this is the raw data obtained from Can
dataFolder = '~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/2019_conc2fluxes/data/';

%names of the files of processed metabolomic data from Leiden
fileName1='20180218 GC-MS Exometabolomics Concentration Calculations_CG.xlsx';
fileName2='20180218 GC-MS Glucose & Lactate Concentration Calculations_CG.xlsx';
fileName3='20180218 LC-MS Exometabolomics Concentration Calculations_CG.xlsx';

%read in the processed metabolomic data
[~, ~, conc1] = xlsread([dataFolder fileName1],'Concentration Values');
[~, ~, SD1] = xlsread([dataFolder fileName1],'Uncertainty');
platform1 = repmat({'GC-MS'}, 1, size(conc1, 2));

[~, ~, conc2] = xlsread([dataFolder fileName2],'Concentration Values');
[~, ~, SD2] = xlsread([dataFolder fileName2],'Uncertainty');
platform2 = repmat({'GC-MS'}, 1, size(conc2, 2));

[~, ~, conc3t] = xlsread([dataFolder fileName3],'Average Concentrations');
[~, ~, SD3t] = xlsread([dataFolder fileName3],'Uncertainty');
platform3 = repmat({'LC-MS'}, 1, size(conc3t, 2));

%remove NaN rows
conc3=conc3t(1:46,:);
SD3=SD3t(1:46,:);

%concatenate the data into one cell array
conc = [conc1, conc2, conc3];
SD = [SD1, SD2, SD3];
platform = [platform1, platform2, platform3];

%check to make sure the sample names are all the same for each row
[mlt,nlt]=size(conc);

ind=find(strcmp('SAMPLE NAME',conc(1,:)));
for m=2:mlt
    if ~(strcmp(conc{m,ind(1)},conc{m,ind(2)}) && strcmp(conc{m,ind(2)},conc{m,ind(3)}))
        warning(['sample names are not consistent for row' num2str(m) ': ' conc{m,ind(1)} ' ' conc{m,ind(2)} ' ' conc{m,ind(3)}])
    end
end

%check to make sure the sample names are all the same for each row
[mlt,nlt]=size(SD);

ind=find(strcmp('SAMPLE NAME',SD(1,:)));
for m=2:mlt
    if ~(strcmp(SD{m,ind(1)},SD{m,ind(2)}) && strcmp(SD{m,ind(2)},SD{m,ind(3)}))
        warning(['sample names are not consistent for row' num2str(m) ': ' SD{m,ind(1)} ' ' SD{m,ind(2)} ' ' SD{m,ind(3)}])
    end
end

%remove redundant columns
conc = [conc(:,1:ind(2)-1), conc(:,ind(2)+1:ind(3)-1), conc(:,ind(3)+1:nlt)];
SD = [SD(:,1:ind(2)-1), SD(:,ind(2)+1:ind(3)-1), SD(:,ind(3)+1:nlt)];
platform = [platform(:,1:ind(2)-1), platform(:,ind(2)+1:ind(3)-1), platform(:,ind(3)+1:nlt)];

leidenMetAbbr = conc(1, 2:length(conc));

%original identifiers
% Key:  A = well number,
%       D = day sample point
%       P = plate
%       e.g. well A1 - Day 13, plate 2, gives D13_P2_A1
% Also, A always corresponds to seeding 400k cells, while B corresponds to seeding 200k cells
%
% 'D13_P3_B1'
% 'D13_P3_B2'
% 'D13_P3_B3'
% 'D13_P4_B1'
% 'D13_P4_B2'
% 'D13_P4_B3'
% 'D19_P3_B1'
% 'D19_P3_B2'
% 'D19_P3_B3'
% 'D19_P4_B1'
% 'D19_P4_B2'
% 'D19_P4_B3'
% 'D23_P4_B1'
% 'D23_P4_B2'
% 'D23_P4_B3'
% 'D13_P3_A1'
% 'D13_P3_A2'
% 'D13_P3_A3'
% 'D13_P4_A1'
% 'D13_P4_A2'
% 'D13_P4_A3'
% 'D19_P3_A1'
% 'D19_P3_A2'
% 'D19_P3_A3'
% 'D19_P4_A1'
% 'D19_P4_A2'
% 'D19_P4_A3'
% 'D23_P4_A1'
% 'D23_P4_A2'
% 'D23_P4_A3'
% 'FRESH MEDIUM 10-12-14 1'
% 'FRESH MEDIUM 10-12-14 2'
% 'FRESH MEDIUM 14-12-14'
% 'MG-MEDIUM A1'
% 'MG-MEDIUM B1'
% 'FRESH MEDIUM D19 1'
% 'FRESH MEDIUM D19 2'
% 'SPENT MEDIUM D19 1'
% 'SPENT MEDIUM D19 2'
% 'P2 A1 13 DAY'
% 'P2 A2 13 DAY'
% 'P2 A3 13 DAY'
% 'P2 B1 13 DAY'
% 'P2 B2 13 DAY'
% 'P2 B3 13 DAY'


%replace identifiers with something consistent and adaptable to Recon
OLD={'P2 A1 13 DAY', 'P2 A2 13 DAY', 'P2 A3 13 DAY', 'P2 B1 13 DAY', 'P2 B2 13 DAY', 'P2 B3 13 DAY',...
    'Pyruvic acid (µM)',...
    'Glycolic acid (µM)',...
    '2-Hydroxybutyric acid (µM)',...
    '3-Hydroxypropionic Acid (µM)',...
    '3-Hydroxybutyric acid (µM)',...
    'Methylmalonic acid (µM)',...
    '3-Hydroxyisovaleric acid (µM)',...
    'Succinic acid (µM)',...
    'Glyceric acid (µM)',...
    'Uracil (µM)',...
    'Fumaric acid (µM)',...
    'Glutaric acid (µM)',...
    'Malic acid (µM)',...
    'Aspartic acid (µM)',...
    'alpha-Ketoglutaric acid (µM)',...
    'Aconitic acid (µM)',...
    'Citric acid (µM)',...
    'Lactic acid (µM)',...
    'Glucose (µM)',...
    'Histidine (µM)',...
    '4-Hydroxyproline (µM)',...
    'Asparagine (µM)',...
    'Taurine (µM)',...
    'Arginine (µM)',...
    'Serine (µM)',...
    'Glutamine (µM)',...
    'Glycine (µM)',...
    'Ethanolamine (µM)',...
    'Aspartic acid (µM)',...
    'Glutamic acid (µM)',...
    'beta-Alanine (µM)',...
    'Threonine (µM)',...
    'Alanine (µM)',...
    'GABA (µM)',...
    'Norepinephrine (µM)',...
    'Cysteine (µM)',...
    'Proline (µM)',...
    'Epinephrine (µM)',...
    'Ornithine (µM)',...
    'Lysine (µM)',...
    'Dopamine (µM)',...
    'Putrescine (µM)',...
    'Tyrosine (µM)',...
    'Methionine (µM)',...
    'Serotonin (µM)',...
    'Valine (µM)',...
    'Isoleucine (µM)',...
    'Leucine (µM)',...
    'Phenylalanine (µM)',...
    'Tryptophan (µM)',...
    'Levodopa (µM)'};
    
NEW={'D13_P2_A1', 'D13_P2_A2', 'D13_P2_A3', 'D13_P2_B1', 'D13_P2_B2','D13_P2_B3',...
    'EX_pyr[e]',...
    'EX_glyclt[e]',...
    'EX_2hb[e]',...
    'EX_3hpp[e]',...
    'EX_bhb[e]',... 
    'EX_HC00900[e]',...
    'EX_CE2028[e]',...
    'EX_succ[e]',...
    'EX_glyc_R[e]',...
    'EX_ura[e]',...
    'EX_fum[e]',...
    'EX_glutar[e]',...
    'EX_mal_L[e]',... 
    'EX_asp_L[e]',... 
    'EX_akg[e]',... 
    'EX_HC00342[e]',...
    'EX_cit[e]',... 
    'EX_lac_L[e]',... 
    'EX_glc_D[e]',...
    'EX_his_L[e]',...
    'EX_4hpro[e]',...
    'EX_asn_L[e]',...
    'EX_taur[e]',... 
    'EX_arg_L[e]',...
    'EX_ser_L[e]',... 
    'EX_gln_L[e]',... 
    'EX_gly[e]',... 
    'EX_etha[e]',...
    'EX_asp_L[e]',... 
    'EX_glu_L[e]',... 
    'EX_ala_B[e]',... 
    'EX_thr_L[e]',...
    'EX_ala_L[e]',... 
    'EX_4abut[e]',... 
    'EX_nrpphr[e]',... 
    'EX_cys_L[e]',...
    'EX_pro_L[e]',... 
    'EX_adrnl[e]',... 
    'EX_orn[e]',... 
    'EX_lys_L[e]',...
    'EX_dopa[e]',... 
    'EX_ptrc[e]',... 
    'EX_tyr_L[e]',... 
    'EX_met_L[e]',...
    'EX_srtn[e]',... 
    'EX_val_L[e]',... 
    'EX_ile_L[e]',... 
    'EX_leu_L[e]',...
    'EX_phe_L[e]',... 
    'EX_trp_L[e]',... 
    'EX_34dhphe[e]'};

[mlt,nlt]=size(conc);
plt=length(OLD);
for m=1:mlt
    for n=1:nlt
        if ischar(conc{m,n})
            for p=1:plt
                conc{m,n}=strrep(conc{m,n},OLD{p},NEW{p});
            end
        end
    end
end

% Every two days, 50% of the spent medium was replaced by fresh medium,
% however samples for exometabolomic data was only obtained on days 13, 19
% and 23.
%
% We assume that the metabolic exchange flux per cell is equal from day 19 to
% 21, and day 21 to 19, and we take into consideration the 50% media change,
% to estimate the spent medium concentration at day 21, and thereby the exchange fluxes.

% smc19 = media concentration at day 19 (KNOWN)
% fmc19 = fresh media concentration at day 19 (50% fresh media, 50% smc19 (KNOWN) )
% smc21 = media concentration at day 21 (UNKNOWN)
% fmc21 = fresh media concentration at day 21 (50% fresh media, 50% smc21(UNKNOWN))
% smc23 = media concentration at day 23 (KNOWN)
% fm = fresh media (KNOWN)

% v = flux rates from fmc19 to smc21 & fmc21 to smc23 (UNKNOWN; assumed to
% be the same in two period times, assuming to do not see
% metabolic phenotypic changes between the population of two time frames)


% fmc19 = 0.5*smc19 + 0.5*fm
% fmc21 = 0.5*smc21 + 0.5*fm
% Solve for smc21
% Assume v = smc21 - fmc19 =  smc23 - fmc21
%            smc21 - 0.5*smc19 - 0.5*fm   =  smc23 - 0.5*smc21 - 0.5*fm
%            smc21 - 0.5*smc19            =  smc23 - 0.5*smc21
%        1.5*smc21 - 0.5*smc19            =  smc23
%        1.5*smc21                        =  smc23 + 0.5*smc19
%            smc21                        =  (2/3)*smc23 + (1/3)*smc19
%
% Calculate v
% v = smc21 - fmc19
% v = (2/3)*smc23 + (1/3)*smc19 - (1/2)*smc19 - (1/2)*fm
% v = (2/3)*smc23 - (1/6)*smc19 - (1/2)*fm

% v =  smc23 - fmc21
% v =  smc23 - (1/2)*smc21 - (1/2)*fm
% v =  smc23 - (1/2)*((2/3)*smc23 + (1/3)*smc19) - (1/2)*fm
% v =  smc23 - (1/3)*smc23 - (1/6)*smc19 - (1/2)*fm
% v =  (2/3)*smc23 - (1/6)*smc19 - (1/2)*fm

% Calculate the standard deviation of v
% vSD = sqrt( (2/3)^2*smc23SD^2 + (1/6)^2*smc19SD^2 + (1/2)^2*fmSD^2 )
% vSD = sqrt( (4/9)*smc23SD^2   + (1/36)*smc19SD^2  + (1/4)*fmSD^2 )

% Homogenous identifiers (not in order)
% 'D13_P2_A1'
% 'D13_P2_A2'
% 'D13_P2_A3'
% 'D13_P2_B1'
% 'D13_P2_B2'
% 'D13_P2_B3'

% 'D13_P3_B1'
% 'D13_P3_B2'
% 'D13_P3_B3'
% 'D13_P4_B1'
% 'D13_P4_B2'
% 'D13_P4_B3'

% 'D13_P3_A1'
% 'D13_P3_A2'
% 'D13_P3_A3'
% 'D13_P4_A1'
% 'D13_P4_A2'
% 'D13_P4_A3'

% 'D19_P3_A1'
% 'D19_P3_A2'
% 'D19_P3_A3'
% 'D19_P3_B1'
% 'D19_P3_B2'
% 'D19_P3_B3'

smc19_400k={...
    'D19_P4_A1';
    'D19_P4_A2';
    'D19_P4_A3'};

smc19_200k={...
    'D19_P4_B1';
    'D19_P4_B2';
    'D19_P4_B3'};

smc19 = [smc19_400k;smc19_200k];

smc23_400k={...
    'D23_P4_A1';
    'D23_P4_A2';
    'D23_P4_A3'};

smc23_200k={...
    'D23_P4_B1';
    'D23_P4_B2';
    'D23_P4_B3'};

smc23 = [smc23_400k;smc23_200k];

% 'FRESH MEDIUM 10-12-14 1'
% 'FRESH MEDIUM 10-12-14 2'
% 'FRESH MEDIUM 14-12-14'

fm={...
    'MG-MEDIUM A1';
    'MG-MEDIUM B1'};

% 'FRESH MEDIUM D19 1'
% 'FRESH MEDIUM D19 2'
% 'SPENT MEDIUM D19 1'
% 'SPENT MEDIUM D19 2'

%calculate the difference in concentration (not normalised)
[~,nlt]=size(conc);
mlt=length(smc19)+1;

dC=cell(mlt,nlt);
dC(1,:)=conc(1,:);
dC{1,1}='CULTURE';
dC{2,1}='P4_A1_400k';
dC{3,1}='P4_A2_400k';
dC{4,1}='P4_A3_400k';
dC{5,1}='P4_B1_200k';
dC{6,1}='P4_B2_200k';
dC{7,1}='P4_B3_200k';

dCSD=cell(mlt,nlt);
dCSD(1,:)=conc(1,:);
dCSD{1,1}='CULTURE';
dCSD{2,1}='P4_A1_400k';
dCSD{3,1}='P4_A2_400k';
dCSD{4,1}='P4_A3_400k';
dCSD{5,1}='P4_B1_200k';
dCSD{6,1}='P4_B2_200k';
dCSD{7,1}='P4_B3_200k';

% v =  (2/3)*smc23 - (1/6)*smc19 - (1/2)*fm
% vSD = sqrt( (4/9)*smc23SD^2   + (1/36)*smc19SD^2  + (1/4)*fmSD^2 )

for n=2:nlt
    for m=2:mlt
        indsmc19  = find(strcmp(smc19{m-1},conc(:,1)));
        indsmc23  = find(strcmp(smc23{m-1},conc(:,1)));
        indfmA1   = find(strcmp('MG-MEDIUM A1',conc(:,1)));
        indfmB1   = find(strcmp('MG-MEDIUM B1',conc(:,1)));
        dC{m,n} = (2/3)*conc{indsmc23,n} - (1/6)*conc{indsmc19,n} - (1/2)*((conc{indfmA1,n}+conc{indfmB1,n})/2);
        dCSD{m,n} = sqrt( (4/9)*(SD{indsmc23,n}^2) + (1/36)*(SD{indsmc19,n}^2) + (1/8)*(SD{indfmA1,n}^2) + (1/8)*(SD{indfmB1,n}^2));
    end
end

%% 1.1  Load cell growth data

% p19 = population at day 19 (KNOWN)
% p21 = population at day 21 (estimated: stimNeuronalGrowthCurve.m)
% p23 = population at day 23 (estimated: stimNeuronalGrowthCurve.m)

% 400k
growthPath = '~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/NeuronalCellGrowth/';
pathfile = 'cellGrowth_400kseeded.mat';
load(strcat(growthPath, pathfile), 'time', 'value');
% obtain data for days of interest
p19_400k = value(strfind(time,19));
%p21_400k = value(strfind(time,21));
%p23_400k = value(strfind(time,23));

% 200k
pathfile = 'cellGrowth_200kseeded.mat';
load(strcat(growthPath, pathfile), 'time', 'value');
clear growthPath pathfile
% obtain data for days of interest
p19_200k = value(strfind(time,19));
%p21_200k = value(strfind(time,21));
%p23_200k = value(strfind(time,23));

% NOTES:
% cell minWeight
% cell maxWeight
% volume (4ml)
% molar (1L)

cellMinW = 0.0000000004396; %grams
cellMaxW = 0.0000000008413; %grams
cellMeanW = (cellMinW + cellMaxW)/2;

interval = 48; % hours
volume = 0.004; % ml (1L Molar)
sign = -1;

%converts the data to uMol/gDW
factorScale = 0.479542184535408; % 0.142701415473448 (old data)

v=dC;
vSD=dCSD;
for n=2:nlt
    for m=2:4
        v{m,n}=v{m,n}*(volume)*(1/(interval*(p19_400k*cellMeanW)))*factorScale;
        vSD{m,n}=vSD{m,n}*(volume)*(1/(interval*(p19_400k*cellMeanW)))*factorScale;
    end
    for m=5:7
        v{m,n}=v{m,n}*(volume)*(1/(interval*(p19_200k*cellMeanW)))*factorScale;
        vSD{m,n}=vSD{m,n}*(volume)*(1/(interval*(p19_200k*cellMeanW)))*factorScale;
    end
end

if 0
    %replace unit
    [mlt,nlt]=size(v);
    for m=1:mlt
        for n=1:nlt
            if ischar(v{m,n})
                v{m,n}=strrep(v{m,n},'(µM)','(µMol/gDW/hr)');
                vSD{m,n}=strrep(vSD{m,n},'(µM)','(µMol/gDW/hr)');
            end
        end
    end
end

%add median exchange rate
v{end+1,1}='Median(µMol/gDW/hr)';
v{end+1,1}='Median absolute deviation';
v{end+1,1}='(Median absolute deviation)/Median';
v{end+1,1}='Mean(µMol/gDW/hr)';
v{end+1,1}='Standard deviation';
v{end+1,1}='Standard error';
v{end+1,1}='lb';
v{end+1,1}='ub';
[mlt,nlt]=size(v);
for n=2:nlt
    v{mlt-7,n}= median([v{5,n},v{6,n},v{7,n}],'omitnan');
    v{mlt-6,n}= mad([v{5,n},v{6,n},v{7,n}],1);
    v{mlt-5,n}= abs(v{mlt-6,n}/v{mlt-7,n});
    v{mlt-4,n}= mean([v{5,n},v{6,n},v{7,n}]);
    v{mlt-3,n}= sqrt(sum([(1/36)*vSD{5,n}^2,(1/36)*vSD{6,n}^2,(1/36)*vSD{7,n}^2]));
    v{mlt-2,n}= abs(v{mlt-3,n}/v{mlt-4,n});
    v{mlt-1,n}= ((cell2mat(v(11, n)) - cell2mat(v(12, n))) * 1); %lb
    v{mlt,n}= ((cell2mat(v(11, n)) + cell2mat(v(12, n))) * 1); %ub
    %vSD{mlt,n}=median([vSD{5,n},vSD{6,n},vSD{7,n}],'omitnan');
    v{mlt+1,n}=leidenMetAbbr{n-1};
    v{mlt+2,n}=platform{n};
end
%save the original identifiers to be able to backtrack in case something
%was misspecified
v{mlt+1,1}='LeidenID';
v{mlt+2,1}='Platform';

clear cellMaxW cellMeanW cellMinW conc conc1 conc2 conc3 conc3t ...
    dataFolder dC dCSD fileName1 fileName2 fileName3 fm ind indfmA1 ...
    indfmB1 indsmc19 indsmc23 interval m mlt n NEW nlt OLD p p19_200k ...
    p19_400k plt SD SD1 SD2 SD3 SD3t sign smc19 smc19_200k smc19_400k ...
    smc23 smc23_200k smc23_400k time value volume vSD

% save(['~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/2019_conc2fluxes/' datestr(now,30) '_can_metabolomics.mat'])