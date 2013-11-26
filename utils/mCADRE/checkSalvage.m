function status=checkSalvage(model)
%when the model is allowed to use PRPP and guanine or hypoxanthine, test if it can
%make gmp or imp.This is the salvage pathway that non-hepatic tissues use
%for purine synthesis.Not useful when the tissue is known to make purines
%de novo

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

model = addSinkReactions(model,{'prpp[c]'},-5,5);

model = changeRxnBounds(model,'EX_gua(e)',-5,'l');
[temp,gmp_dm]=addDemandReaction(model,'gmp[c]');
temp=changeObjective(temp,gmp_dm,1);
sol=optimizeCbModel(temp);
status_gmp=sol.f>1e-6;

model = changeRxnBounds(model,'EX_gua(e)',0,'l');
model = changeRxnBounds(model,'EX_hxan(e)',-5,'l');
[temp,imp_dm]=addDemandReaction(model,'imp[c]');
temp=changeObjective(temp,imp_dm,1);
sol=optimizeCbModel(temp);
status_imp=sol.f>1e-6;

status=status_gmp&&status_imp;