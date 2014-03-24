function [modelOut, LBorUBdiff] = removeEnzymeIrrevs(model)
% We use the convention that only specified rev rxns
% are actually reversible, meaning that we must relax
% constraints on lb/ub to match this.

[selExc, selUpt] = findExcRxns(model);
modelOut = model;
geneRxns = boolean(sum(model.rxnGeneMat')');
lbRxns = find(model.rev & (model.lb >= 0) & ~selExc & geneRxns);
ubRxns = find(model.rev & (model.ub <= 0) & ~selExc & geneRxns);
modelOut.lb(lbRxns) = -model.ub(lbRxns);
modelOut.ub(ubRxns) = -model.lb(ubRxns);

% say how many bounds have changed.
LBchanged = modelOut.lb ~= model.lb;
UBchanged = modelOut.ub ~= model.ub;
LBorUBdiff = LBchanged | UBchanged;
if isfield(model, 'description')
    disp([model.description ' removed Irrevs:']);
end
disp(sprintf('Lower bounds changed: %d', sum(LBchanged))); 
disp(sprintf('Upper bounds changed: %d', sum(UBchanged))); 
disp(sprintf('Total bounds changed: %d', sum(LBorUBdiff))); 