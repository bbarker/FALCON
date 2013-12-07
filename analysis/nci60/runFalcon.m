function [v_solirrev v_solrev corral] = ...
    runFalcon(model, expressionFile, rc, EXPCON, FDEBUG)
%INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   lb           Lower bounds
%   ub           Upper bounds
%   grRules      boolean rules in terms of genes (see next field)
%   genes        Entrez gene ids (without isoform id)
%
% expressionFile    Tab-delimited file with a head for the columns:
%                   gene (entrez gene id), mean (expression value,
%                   and standard deviation (of expression).
%
% rc                regularization constant on fluxes
%
% EXPCON            Flag to specify whether to use expression constraints.
%
%OUTPUT
% v_solirrev    irreversible flux vector
%
% v_solrev      reversible flux vector
%
%
% Yiping Wang    09/08/13
% Brandon Barker 09/15/2013  Now calls convertIrrevFluxDistribution

if ~exist('FDEBUG', 'var')
    FDEBUG = false;
end

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

[rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
    computeMinDisj(modelIrrev, expressionFile);

[v_solirrev corral]= falcon(modelIrrev, rxn_exp, rxn_exp_sd, ...
                    rxn_rule_group, rc, 0, EXPCON, FDEBUG);
for i=1:length(modelIrrev.rxns)
  %disp([modelIrrev.rxns{i} ' ' num2str(v_solirrev(i))]);
end
v_solrev = convertIrrevFluxDistribution(v_solirrev, matchRev);

curatedRxnNames={'HEX1','PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','LDH_L','G6PDH2c','G6PDH2r','PGL','PGLc','GND','GNDc','RPI','RPE','RPEc','TALA','TKT1','TKT2','PRPPS','GLUPRT','PRAGSr','GARFT','PRFGS','r0666','AIRCr','PRASCS','ADSL2','AICART','IMPC'};
glycolysisRxnNames=model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis'));
pentosePhosphateRxnNames=model.rxns(strcmp(model.subSystems,'Pentose phosphate pathway'));
purineSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Purine synthesis'));
citricAcidCycleRxnNames=model.rxns(strcmp(model.subSystems,'Citric acid cycle'));
oxidativePhosphorylationRxnNames=model.rxns(strcmp(model.subSystems,'Oxidative phosphorylation'));
curatedRxnNames=[glycolysisRxnNames; pentosePhosphateRxnNames; purineSynthesisRxnNames; oxidativePhosphorylationRxnNames; citricAcidCycleRxnNames];
curatedRxnFluxes=zeros(1,length(curatedRxnNames));
outputFI=fopen('../rxnsToShow.txt','w');
for i=1:length(model.rxns)
  curatedRxnIdx=find(strcmp(model.rxns{i},curatedRxnNames));
  if(sum(strcmp(model.rxns{i},curatedRxnNames))~=0)
    disp([model.rxns{i} ' ' num2str(v_solrev(i))]);
  end
  if(sum(strcmp(model.rxns{i},curatedRxnNames))==0 && abs(v_solrev(i))>0 && sum(strfind(model.rxns{i},'TR_'))==0)
    %disp([model.rxns{i} ' ' num2str(v_solrev(i))]);
    %fprintf(outputFI,'%s\n',model.rxns{i});
  end
  curatedRxnFluxes(curatedRxnIdx)=v_solrev(i);
  if(sum(strfind(model.rxns{i},'TR_'))~=0)
    disp([model.rxns{i} ' ' num2str(v_solrev(i))]);
  end
end

%makeBarGraph(1:length(curatedRxnNames),curatedRxnFluxes,1:length(curatedRxnNames),curatedRxnNames,[-1000 1000],'curatedRxnFluxes','../curatedRxnFluxes.png');