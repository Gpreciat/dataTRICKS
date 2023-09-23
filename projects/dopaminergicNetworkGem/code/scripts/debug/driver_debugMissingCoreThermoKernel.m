% --- thermoKernel END ----
% 5476 reactions removed by createTissueSpecificModel.
% 26 core reactions removed by createTissueSpecificModel.
%     <strong>Forward_Reaction</strong>                           <strong>Name</strong>                           <strong>lb</strong>      <strong>ub</strong>                                                                                                         <strong>equation</strong>                                                                                                    
%     <strong>________________</strong>    <strong>__________________________________________________</strong>    <strong>__</strong>    <strong>______</strong>    <strong>_______________________________________________________________________________________________________________________________________________________________________________________________________________</strong>
% 
%     {'ALOX5'       }    {'Arachidonate 5-Lipoxygenase'                   }    0     100000    {'o2[c] + arachd[c]  -> 5HPET[c] '                                                                                                                                                                            }
%     {'AMPDA'       }    {'Adenosine Monophosphate Deaminase'             }    0     100000    {'h2o[c] + h[c] + amp[c]  -> nh4[c] + imp[c] + dummy_Met_271 '                                                                                                                                                }
%     {'DAGK_hs'     }    {'Diacylglycerol Phosphate Kinase (Homo Sapiens)'}    0     100000    {'atp[c] + dag_hs[c]  -> h[c] + adp[c] + pa_hs[c] + dummy_Met_1606 + dummy_Met_1607 + dummy_Met_1608 + dummy_Met_160851 + dummy_Met_1609 + dummy_Met_8525 + dummy_Met_8526 + dummy_Met_8527 + dummy_Met_9162 '}
%     {'EX_co[e]'    }    {'Exchange of Carbon Monoxide '                  }    0      10000    {'co[e]  -> '                                                                                                                                                                                                 }
%     {'EX_lneldc[e]'}    {'Exchange of Linoelaidic Acid '                 }    0     100000    {'lneldc[e]  -> '                                                                                                                                                                                             }
%     {'FPGS6'       }    {'Folylpolyglutamate Synthetase (Dhf)'           }    0     100000    {'atp[c] + glu_L[c] + 6dhf[c]  -> h[c] + adp[c] + pi[c] + 7dhf[c] + dummy_Met_2356 '                                                                                                                          }
%     {'FPGS8'       }    {'Folylpolyglutamate Synthetase (10Fthf)'        }    0     100000    {'10fthf5glu[c] + atp[c] + glu_L[c]  -> 10fthf6glu[c] + h[c] + adp[c] + pi[c] + dummy_Met_2356 '                                                                                                              }
%     {'FPGS9'       }    {'Folylpolyglutamate Synthetase (10Fthf)'        }    0     100000    {'10fthf6glu[c] + atp[c] + glu_L[c]  -> 10fthf7glu[c] + h[c] + adp[c] + pi[c] + dummy_Met_2356 '                                                                                                              }
%     {'NMNATn'      }    {'Nicotinamide-Nucleotide Adenylyltransferase'   }    0     100000    {'atp[n] + h[n] + nmn[n]  -> ppi[n] + nad[n] + dummy_Met_64802 '                                                                                                                                              }
%     {'PLA2_2'      }    {'Phospholipase A2'                              }    0     100000    {'h2o[c] + pchol_hs[c]  -> h[c] + Rtotal2[c] + lpchol_hs[c] + dummy_Met_5321 + dummy_Met_81579 + dummy_Met_8398 '                                                                                             }
% 
%     <strong>Reversible_Reaction</strong>                         <strong>Name</strong>                                  <strong>lb</strong>               <strong>ub</strong>                                    <strong>equation</strong>                               
%     <strong>___________________</strong>    <strong>______________________________________________</strong>    <strong>____________________</strong>    <strong>______</strong>    <strong>_____________________________________________________________________</strong>
% 
%      {'EX_ascb_L[e]' }     {'Exchange of L-Ascorbate '                  }       -33.4428226630558    100000    {'ascb_L[e]  <=> '                                                  }
%      {'EX_lipoate[e]'}     {'Exchange of Lipoate '                      }     -0.0394044564234044    100000    {'lipoate[e]  <=> '                                                 }
%      {'EX_lnlc[e]'   }     {'Exchange of Linoleic Acid (All Cis C18:2) '}     -0.0115961686929536    100000    {'lnlc[e]  <=> '                                                    }
%      {'EX_ncam[e]'   }     {'Exchange of Nicotinamide '                 }       -3.81469600757283    100000    {'ncam[e]  <=> '                                                    }
%      {'EX_prgstrn[e]'}     {'Exchange of Progesterone '                 }    -0.00161329356250689    100000    {'prgstrn[e]  <=> '                                                 }
%      {'RE1530M'      }     {'Thymidine Kinase'                          }                 -100000    100000    {'dgtp[m] + duri[m]  <=> h[m] + dgdp[m] + dump[m] + dummy_Met_7084 '}
%      {'RE1918C'      }     {'RE1918C'                                   }                 -100000    100000    {'dopa[c] + acald[c]  <=> h2o[c] + C09642[c] '                      }
%      {'EX_pydxn[e]'  }     {'Exchange of Pyridoxine'                    }       -2.26640035357874    100000    {'pydxn[e]  <=> '                                                   }
% 
% 
% 3458 metabolites removed by createTissueSpecificModel.
% 4 core metabolites removed by createTissueSpecificModel.
%     {'phyQ[c]' }
%     {'retfa[c]'}
%     {'sprm[c]' }
%     {'thmpp[c]'}


load('~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/iDopaNeuro1/19.debug_prior_to_create_tissue_specific_model.mat')

%%
if 1
    %remove dummy reactions and metabolites
    model = destroyDummyModel(model);
    % remove coupling constraints
    model = removeCOBRAConstraints(model, model.ctrs);
end

%%
coreMetAbbr = {'phyQ[c]','retfa[c]','sprm[c]','thmpp[c]'};
coreRxnAbbr = {'ALOX5','AMPDA','DAGK_hs','EX_co[e]','EX_lneldc[e]','FPGS6','FPGS8','FPGS9','NMNATn','PLA2_2','EX_ascb_L[e]','EX_lipoate[e]','EX_lnlc[e]',...
    'EX_ncam[e]','EX_prgstrn[e]','RE1530M','RE1918C','RE1918C','EX_pydxn[e]'}; 

%%
setdiff(coreMetAbbr,model.mets)
setdiff(coreRxnAbbr,model.rxns)

%%
paramThermoFluxConsistency.formulation = 'pqzw';
paramThermoFluxConsistency.epsilon = param.thermoFluxEpsilon;
paramThermoFluxConsistency.printLevel = param.printLevel-1;
paramThermoFluxConsistency.nMax = 50;
paramThermoFluxConsistency.relaxBounds=0;
paramThermoFluxConsistency.acceptRepairedFlux=1;
paramThermoFluxConsistency.iterationMethod = 'random';
[thermoFluxConsistentMetBool, thermoFluxConsistentRxnBool, model, thermoConsistModel]...
    = findThermoConsistentFluxSubset(model, paramThermoFluxConsistency);

%%
model = thermoConsistModel;

%%
setdiff(coreMetAbbr,thermoConsistModel.mets)
setdiff(coreRxnAbbr,thermoConsistModel.rxns)

%%
load('~/work/sbgCloud/programReconstruction/projects/exoMetDN/results/thermoKernelDebug/model_to_debug')

%%
tissueModelOptions.metWeights = zeros(length(model.mets), 1);
tissueModelOptions.metWeights(ismember(model.mets,coreMetAbbr))=-1;

tissueModelOptions.rxnWeights = zeros(length(model.rxns), 1);
tissueModelOptions.rxnWeights(ismember(model.rxns,coreRxnAbbr))=-1;

tissueModelOptions.solver = 'thermoKernel';
tissueModelOptions.printLevel = param.printLevel-1;
tissueModelOptions.formulation = 'pqzwrs';
tissueModelOptions.nMax = 20;
tissueModelOptions.relaxBounds = 0;
tissueModelOptions.acceptRepairedFlux = 1;
%tissueModelOptions.iterationMethod = 'greedyRandom';
tissueModelOptions.iterationMethod = 'random';
%tissueModelOptions.iterationMethod = 'greedy';
%tissueModelOptions.iterationMethod = 'greedyRandom2';
tissueModelOptions.normalizeZeroNormWeights = 0;
tissueModelOptions.epsilon = param.thermoFluxEpsilon;

modelTemp = createTissueSpecificModel(model, tissueModelOptions, 0);

%%
% Look for removed reactions
bool = ~ismember(model.rxns, modelTemp.rxns);
rxnsNotInModel = model.rxns(bool);
if ~isempty(rxnsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n', length(rxnsNotInModel), ' reactions removed by createTissueSpecificModel.')
    end
    if param.printLevel > 2
        printConstraints(model, -inf, inf, bool)
    end
end
% Look for core reactions removed
coreRxnsNotInModel = setdiff(coreRxnAbbr,modelTemp.rxns);
if ~isempty(coreRxnsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n',length(coreRxnsNotInModel), ' core reactions removed by createTissueSpecificModel.')
        if param.printLevel > 1
            printConstraints(model, -inf, inf, ismember(model.rxns,coreRxnsNotInModel))
        end
    end
end
% Look for removed metabolites
bool = ~ismember(model.mets, modelTemp.mets);
metsNotInModel = model.mets(bool);
if ~isempty(metsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n', length(metsNotInModel), ' metabolites removed by createTissueSpecificModel.')
    end
end
% Look for core metabolites removed
coreMetsNotInModel = setdiff(coreMetAbbr,modelTemp.mets);
if ~isempty(coreMetsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n', length(coreMetsNotInModel), ' core metabolites removed by createTissueSpecificModel.')
        if param.printLevel
            disp(coreMetsNotInModel)
        end
    end
end
return

%%
coreRxnAbbr = {'ALOX5','AMPDA','DAGK_hs','EX_co[e]','EX_lneldc[e]','FPGS6','FPGS8','FPGS9','NMNATn','PLA2_2','EX_lipoate[e]','EX_lnlc[e]',...
    'EX_ncam[e]','EX_prgstrn[e]','RE1530M','RE1918C','RE1918C','EX_pydxn[e]'}; 

coreMetAbbr = {'phyQ[c]','retfa[c]','thmpp[c]'};


tissueModelOptions.metWeights = zeros(length(model.mets), 1);
tissueModelOptions.metWeights(ismember(model.mets,coreMetAbbr))=-1;

tissueModelOptions.rxnWeights = zeros(length(model.rxns), 1);
tissueModelOptions.rxnWeights(ismember(model.rxns,coreRxnAbbr))=-1;

modelTemp = createTissueSpecificModel(model, tissueModelOptions, 0);

% Look for removed reactions
bool = ~ismember(model.rxns, modelTemp.rxns);
rxnsNotInModel = model.rxns(bool);
if ~isempty(rxnsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n', length(rxnsNotInModel), ' reactions removed by createTissueSpecificModel.')
    end
    if param.printLevel > 2
        printConstraints(model, -inf, inf, bool)
    end
end
% Look for core reactions removed
coreRxnsNotInModel = setdiff(coreRxnAbbr,modelTemp.rxns);
if ~isempty(coreRxnsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n',length(coreRxnsNotInModel), ' core reactions removed by createTissueSpecificModel.')
        if param.printLevel > 1
            printConstraints(model, -inf, inf, ismember(model.rxns,coreRxnsNotInModel))
        end
    end
end
% Look for removed metabolites
bool = ~ismember(model.mets, modelTemp.mets);
metsNotInModel = model.mets(bool);
if ~isempty(metsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n', length(metsNotInModel), ' metabolites removed by createTissueSpecificModel.')
    end
end
% Look for core metabolites removed
coreMetsNotInModel = setdiff(coreMetAbbr,modelTemp.mets);
if ~isempty(coreMetsNotInModel)
    if param.printLevel > 0
        fprintf('%u%s\n', length(coreMetsNotInModel), ' core metabolites removed by createTissueSpecificModel.')
        if param.printLevel
            disp(coreMetsNotInModel)
        end
    end
end

%%
paramThermoFluxConsistency.formulation = 'pqzw';
paramThermoFluxConsistency.epsilon = param.thermoFluxEpsilon;
paramThermoFluxConsistency.printLevel = param.printLevel-1;
paramThermoFluxConsistency.nMax = 40;
paramThermoFluxConsistency.relaxBounds=0;
paramThermoFluxConsistency.acceptRepairedFlux=1;
paramThermoFluxConsistency.iterationMethod = 'random';
[thermoFluxConsistentMetBool, thermoFluxConsistentRxnBool, model, thermoConsistModel]...
    = findThermoConsistentFluxSubset(model, paramThermoFluxConsistency);

%%
load('debug_thermoKernel');
%%
options.metWeights(:)=0;
%options.rxnWeights(options.rxnWeights>0)=0;
[tissueModel, ~, ~] = thermoKernel(model, options.activeInactiveRxn, options.rxnWeights, options.presentAbsentMet, options.metWeights, options);