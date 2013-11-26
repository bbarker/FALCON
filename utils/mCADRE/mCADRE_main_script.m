load indGeneRxn %indices of gene-associated reactions in the human Recon 1
load reactionScores %expression-based evidence scores for gene-associated reactions. Should be in the same order as indGeneRxn
load model %human recon1 with the original gene id converted to numerical entrez gene id
load rxn_map_matrix %the adjacency matrix of reactions connected by common metabolites. Diagonal elements are 1
load confidenceScores %confidence scores for all reactions in recon1
load precursorMets %metabolites that should be produced from glucose by the tissue-specific model
load HR1_CbModel %human recon 1 with the original gene id format
initCobraToolbox

model.genes=cellstr(num2str(model.genes));%convert numerical entrez gene id to cell array of strings

numrxns=numel(model.rxns);
reactionScoresAll=zeros(numrxns,1);
reactionScoresAll(indGeneRxn,1)=reactionScores(:,GXX);
reactionScoresAll(isnan(reactionScoresAll))=0;
%assign expression-based evidence to gene-associated reactions

core=model.rxns(reactionScoresAll>=0.5);
%core reactions defined as reactions present in at least 50% of samples
%based on data. Can be a different threshold.

inactive_g=checkModelConsistency(model,[],[],1);
active_g=setdiff(model.rxns,inactive_g);
core_active_g=intersect(core,active_g);
%identify core reactions that can also carry flux
 
ind2active=ismember(model.rxns,active_g);
active_g=model.rxns(ind2active);
rxn_map_matrix=rxn_map_matrix-eye(size(rxn_map_matrix));
connect2core=rxn_map_matrix(ind2active,ind2active);
connect2core=connect2core./repmat(sum(connect2core,2),1,size(connect2core,2));
connect2coreScore=connect2core*reactionScoresAll(ind2active,1);
connect2coreScore(isnan(connect2coreScore))=0;
%calculate the connection-based evidence

non_core=setdiff(model.rxns,union(core,setdiff(model.rxns,active_g)));
[c,ind]=ismember(non_core,active_g);
connect=connect2coreScore(ind(c));
[c,ind]=ismember(non_core,model.rxns);
exp_score=reactionScoresAll(ind(c),1);
confidence=confidenceScores(ind(c),1);
generxns=model.rxns(indGeneRxn);
ind2GeneNoExp=ismember(non_core,generxns);
exp_score(ind2GeneNoExp,1)=exp_score(ind2GeneNoExp,1)-1e-6;
%do this to distinguish gene associated reaction NOT expressed (0 score)
%and non-gene associated reaction (0 score). In removal,gene
%associated reaction NOT expressed to be removed first. 
[c,indNonCore]=sortrows([connect,exp_score,confidence],[2,1,3]);
%sort non_core rxns first by connectivity to core, then by expression, then
%by confidenceScores in Recon1.
rxns2remove=non_core(indNonCore);
%the above section map connection-based,confidence score-based
%expression-based scores to non-core reactions and rank them

count=0;
tissue_model=removeRxns(model,setdiff(model.rxns,active_g));
%the initial tissue_model will be a model with all the 2469 flux-carrying
%reactions in recon1

precursorMets=[precursorMets;nonEssentialAA;nucleotide;lipid];
if checkModelFunction(tissue_model,precursorMets)
display('genericModel passed precursor metabolites test!');
%check if the generic model can pass the metabolite production test
while ~isempty(rxns2remove)
    rxn2remove=rxns2remove(1);
    inactive=checkModelConsistency(tissue_model,rxn2remove,core_active_g,1);%reactions that are inactivated by the removal of the current non-core reaction.
    if exp_score(ismember(non_core,rxn2remove))==-1e-6&&... %if the non-core reaction is not present in ANY tissue samples
    numel(intersect(inactive,core_active_g))/numel(intersect(inactive,non_core))<=(1/3)&&... %if the number of non-core reactions being removed exceeds the number of core reactions being removed exceeds 3 to 1
    checkModelFunction(removeRxns(tissue_model,inactive),precursorMets)&&checkSalvage(removeRxns(tissue_model,inactive)) %and if removing these reactions won't affect tissue model's core functinalities
    %checkSalvage tests the ability of the tissue model to use the salvage pathway to make purines, as de novo purine synthesis primarily occurs in the liver.  Not useful when the tissue is known to make purines de novo  
        tissue_model=removeRxns(tissue_model,inactive);
        rxns2remove=rxns2remove(~ismember(rxns2remove,inactive),1);
        display(['Number of rxns removed in this iter: ',num2str(numel(inactive))])    
        count=count+numel(inactive);   
        display(['Number of rxns removed so far: ',num2str(count)])
        core_inactivated=intersect(core,inactive);        
        core=setdiff(core,inactive);
        core_active_g=setdiff(core_active_g,inactive);
        display(['Number of core rxns removed: ',num2str(numel(core_inactivated))])   
        display(['Number of rxns to remove: ',num2str(numel(rxns2remove))])  
        continue;
    end
    %Reactions not present in any of the tissue samples based on expression
    %data form a high confidence negative set. mCADRE tries to remove them
    %even if some core reactions have to be removed, if 3 conditions are
    %met, as commented above
    partofcore=intersect(inactive,core_active_g);
    inactive=setdiff(inactive,core);
    temp_model=removeRxns(tissue_model,inactive);
    if numel(partofcore)==0&&checkModelFunction(temp_model,precursorMets)&&checkSalvage(temp_model)
        tissue_model=temp_model;
        rxns2remove=rxns2remove(~ismember(rxns2remove,inactive),1);%changed on May-09-2012
        display(['Number of rxns removed in this iter: ',num2str(numel(inactive))])    
        count=count+numel(inactive);
    else
        rxns2remove(1)=[];
        display(['Number of rxns removed in this iter: 0'])
    end
    display(['Number of rxns removed so far: ',num2str(count)])
    display(['Number of rxns to remove: ',num2str(numel(rxns2remove))])
    %For all other non-core reactions, flux through the core reactions, and
    %core model functionality must be retained.
end

rxnsInTissue=tissue_model.rxns;
tissue_model_originalID=removeRxns(originalModel,setdiff(originalModel.rxns,rxnsInTissue));
%tissue_model_originalID retains the original gene id format in recon1

tissue_model=removeNonUsedGenes(tissue_model);
tissue_model_originalID=removeNonUsedGenes(tissue_model_originalID);
%removeNonUsedGenes removes genes no longer used in gene-reaction rules
%after the corresponding reactions are removed.
save TissueReconResultsGXX tissue_model tissue_model_originalID
else
    display('genericModel failed to pass precursor metabolites test!');
end