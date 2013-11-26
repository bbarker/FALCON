function [precursorMetsStatus,time] = checkModelFunction(model,precursorMets)

t0 = clock;

indExRxns = findExRxns(model);
exRxns = model.rxns(indExRxns);

carbonMets = regexp(model.metFormulas,'C');
hydrogenMets = regexp(model.metFormulas,'H');
is_organic = ~cellfun('isempty',carbonMets)&~cellfun('isempty',hydrogenMets);
organicMets = model.mets(is_organic);

organicRxns = findRxnsFromMets(model,organicMets);
organicExRxns = intersect(organicRxns,exRxns);
organicExRxns = [organicExRxns;...
    'EX_Rtotal(e)';'EX_Rtotal2(e)';'EX_Rtotal3(e)';'EX_Tyr_ggn(e)'];
model = changeRxnBounds(model,organicExRxns,0,'l');

model = changeRxnBounds(model,'EX_glc(e)',-5,'l');
model = changeRxnBounds(model,'EX_co2(e)',-1000,'l');
 
coaMets= {'accoa[m]';'succoa[m]';'pmtcoa[c]'};
[model,addedRxns] = addDemandReaction(model,setdiff(precursorMets,coaMets));
model = addReaction(model,'DM_accoa(m)',{'accoa[m]';'coa[m]'},[-1;1],0,0,1000,0);
model = addReaction(model,'DM_succoa(m)',{'succoa[m]';'coa[m]'},[-1;1],0,0,1000,0);
model = addReaction(model,'DM_pmtcoa(c)',{'pmtcoa[c]';'coa[c]'},[-1;1],0,0,1000,0);
precursorRxns = [addedRxns,'DM_accoa(m)','DM_succoa(m)','DM_pmtcoa(c)']';

rxnList = precursorRxns;
inactivePrecursors = [];
while numel(rxnList)
    numRxnList = numel(rxnList);
    model = changeObjective(model,rxnList);
    
    % Maximize all
    FBAsolution = optimizeCbModel(model,'max');
    optMax = FBAsolution.x;
    if isempty(optMax)
        inactivePrecursors=1;
        break;
    end
    precursorFlux = optMax(ismember(model.rxns,precursorRxns));
    activePrecursors = precursorRxns(abs(precursorFlux) >= 1e-6);
    rxnList = setdiff(rxnList,activePrecursors);
    
    numRemoved = numRxnList - numel(rxnList);
    
    if ~numRemoved
        randInd = randperm(numel(rxnList));
        i = rxnList(randInd(1));
        model = changeObjective(model,i);
        
        % Maximize reaction i
        FBAsolution = optimizeCbModel(model,'max');
        optMax = FBAsolution.f;
        if abs(optMax) < 1e-6
            inactivePrecursors = union(inactivePrecursors,i);
            break;
        end
        
        rxnList = setdiff(rxnList,i);
    end
end

precursorMetsStatus = ~numel(inactivePrecursors);
time = etime(clock,t0);
