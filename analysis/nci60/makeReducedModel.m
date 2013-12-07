function reducedModel=makeReducedModel(model)
%reducedModel=extractSubNetwork(model,{'HEX1','PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK',model.rxns{strcmp(model.subSystems,'Exchange/demand reaction')},model.rxns{strcmp(model.subSystems,'Transport, extracellular')},model.rxns{strcmp(model.subSystems,'Transport, mitochondrial')},model.rxns{strcmp(model.subSystems,'Citric acid cycle')},model.rxns{strcmp(model.subSystems,'Oxidative phosphorylation')}});
%reducedModel=extractSubNetwork(model,{model.rxns{strcmp(model.subSystems,'Glycolysis/gluconeogenesis')},model.rxns{strcmp(model.subSystems,'Exchange/demand reaction')},model.rxns{strcmp(model.subSystems,'Transport, extracellular')},model.rxns{strcmp(model.subSystems,'Transport, mitochondrial')},model.rxns{strcmp(model.subSystems,'Citric acid cycle')},model.rxns{strcmp(model.subSystems,'Oxidative phosphorylation')}});

curatedRxnNames={'HEX1','PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','LDH_L','G6PDH2c','G6PDH2r','PGL','PGLc','GND','GNDc','RPI','RPE','RPEc','TALA','TKT1','TKT2','PRPPS','GLUPRT','PRAGSr','GARFT','PRFGS','r0666','AIRCr','PRASCS','ADSL2','AICART','IMPC','PYRt2m'};

simpleGlycolysisRxnNames=curatedRxnNames;
glycolysisRxnNames=model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis'));
pentosePhosphateRxnNames=model.rxns(strcmp(model.subSystems,'Pentose phosphate pathway'));
purineSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Purine synthesis'));
citricAcidCycleRxnNames=model.rxns(strcmp(model.subSystems,'Citric acid cycle'));
oxidativePhosphorylationRxnNames=model.rxns(strcmp(model.subSystems,'Oxidative phosphorylation'));
exchangeDemandRxnNames=model.rxns(strcmp(model.subSystems,'Exchange/demand reactions'));
transportRxnNames=model.rxns(strcmp(model.subSystems,'Transport, mitochondrial'));
simpleGlycolysisRxnNames=[glycolysisRxnNames; pentosePhosphateRxnNames; purineSynthesisRxnNames; transportRxnNames; oxidativePhosphorylationRxnNames; citricAcidCycleRxnNames];%; exchangeDemandRxnNames];%; simpleGlycolysisRxnNames']
disp(transportRxnNames{257});
selRxns = ismember(model.rxns,simpleGlycolysisRxnNames);
subS = model.S(:,selRxns);
if (nargin < 3)
    selMets = ~all(subS == 0,2);
else
    selMets = ismember(model.mets,metNames);
end

subS = subS(selMets,:);

subModel.S = subS;
subModel.rxns = model.rxns(selRxns);
subModel.mets = model.mets(selMets);
if (isfield(model,'b'))
    subModel.b = model.b(selMets);
end
if (isfield(model,'metNames'))
    subModel.metNames = model.metNames(selMets);
end
if (isfield(model,'metFormulas'))
    subModel.metFormulas = model.metFormulas(selMets);
end
if (isfield(model,'description'))
    subModel.description = model.description;
end
if (isfield(model,'rev'))
    subModel.rev = model.rev(selRxns);
end
if (isfield(model,'lb'))
    subModel.lb = model.lb(selRxns);
end 
if (isfield(model,'ub'))
    subModel.ub = model.ub(selRxns);
end
if (isfield(model,'c'))
    subModel.c = model.c(selRxns);
end
if (isfield(model,'genes'))
   newRxnGeneMat = model.rxnGeneMat(selRxns,:); 
   selGenes = sum(newRxnGeneMat)' > 0;
   subModel.rxnGeneMat = newRxnGeneMat(:,selGenes);
   subModel.genes = model.genes(selGenes);
   subModel.grRules = model.grRules(selRxns);
   subModel.rules = model.rules(selRxns);
end
if (isfield(model,'geneNames'))
    subModel.geneNameRules = model.geneNameRules(selRxns);
    subModel.geneNames = model.geneNames(selGenes);
end
if (isfield(model,'subSystems'))
    subModel.subSystems = model.subSystems(selRxns);
end
if(isfield(model,'rxnNames'))
    subModel.rxnNames=model.rxnNames(selRxns);
end
if isfield(model,'reversibleModel')
    subModel.reversibleModel=model.reversibleModel;
end

reducedModel=subModel;
transportMets={'h2o[c]','h[c]','pi[c]','lac_L[c]','glc_D[c]','co2[c]','gln_L[c]','glu_L[c]','gly[c]','asp_L[c]','nadp[c]','nadph[c]','atp[c]','adp[c]','amp[c]','ppi[c]','10fthf[c]','thf[c]','fum[c]','imp[c]'};
transportRxns={}; 

%pentosePhosphateRxnNames={'PGL','PGLc','GND','GNDc','RPI','RPE','RPEc','TALA','TKT1','TKT2','PRPPS','GLUPRT','PRAGSr','GARFT','PRFGS','r0666','AIRCr','PRASCS','ADSL2','AICART','IMPC','PYRt2m'};

for i=1:length(reducedModel.rxns)
    if(sum(strcmp(reducedModel.rxns{i},pentosePhosphateRxnNames))~=0)
        %reducedModel.rev(i)=1;
	%reducedModel.lb(i)=-1000;
	%reducedModel.ub(i)=1000;
    end
    %reducedModel.rev(i)=1;
    %reducedModel.lb(i)=-1000;
    %reducedModel.ub(i)=1000;
    if(strcmp(reducedModel.rxns{i},'PYK'))
        reducedModel.lb(i)=1;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=0;
    end
    if(strcmp(reducedModel.rxns{i},'IMPC'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=-1;
	reducedModel.rev(i)=0;
    end
end

for i=1:length(reducedModel.mets)
    %internalMet=[reducedModel.mets{i} '[c]'];
    %externalMet=[reducedModel.mets{i} '[e]'];
    %internalIdx=strcmp(reducedModel.mets,internalMet);
    %externalIdx=strcmp(reducedModel.mets,externalMet);
    internalIdx=i;
    noConsumingReaction=sum(reducedModel.S(internalIdx,:)<0)==0;
    noProducingReaction=sum(reducedModel.S(internalIdx,:)>0)==0;
    %noConsumingReaction
    %if((noProducingReaction && noConsumingReaction) || sum(strcmp([reducedModel.mets{i}],transportMets))~=0)
    if(sum(strcmp([reducedModel.mets{i}],transportMets))~=0)
    reducedModel.rxns{end+1}=['TR_' reducedModel.mets{i}];
    reducedModel.S(internalIdx,end+1)=1;
    %reducedModel.S(externalIdx,end)=-1;
    reducedModel.rev(end+1)=1;
    reducedModel.lb(end+1)=-10000;
    reducedModel.ub(end+1)=10000;
    reducedModel.c(end+1)=0;
    reducedModel.rxnGeneMat(end+1,:)=zeros(1,size(reducedModel.rxnGeneMat,2));
    reducedModel.grRules{end+1}='';
    reducedModel.subSystems{end+1}='';
    reducedModel.rxnNames{end+1}='';
    reducedModel.rules{end+1}='';
    end
end



if(0)
reducedModel.rxns{end+1}='D_nadh';
internalMet='nadh[c]';
%externalMet=[transportMets{i} '[e]'];
internalIdx=strcmp(reducedModel.mets,internalMet);
%externalIdx=strcmp(reducedModel.mets,externalMet);
reducedModel.S(internalIdx,end+1)=-1;
%reducedModel.S(externalIdx,end)=-1;
reducedModel.rev(end+1)=1;
reducedModel.lb(end+1)=-1000;
reducedModel.ub(end+1)=1000;
reducedModel.c(end+1)=0;
reducedModel.rxnGeneMat(end+1,:)=zeros(1,size(reducedModel.rxnGeneMat,2));
reducedModel.grRules{end+1}='';
reducedModel.subSystems{end+1}='';
reducedModel.rxnNames{end+1}='';
reducedModel.rules{end+1}='';
end