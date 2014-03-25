function [modelOut, revToIrrev] = setEnzymeIrrevs(model)
% This script does the oppposite of removeEnzymeIrrevs, and
% sets the rev field to match the lb and ub constraints.

% For now this just does unirectional -> set irrev, not
% bidirectional -> set rev.

[selExc, selUpt] = findExcRxns(model);
modelOut = model;
geneRxns = boolean(sum(model.rxnGeneMat')');
lbRxns = find(model.rev & (model.lb >= 0) & ~selExc & geneRxns);
ubRxns = find(model.rev & (model.ub <= 0) & ~selExc & geneRxns);
%modelOut.lb(lbRxns) = -model.ub(lbRxns); % in removeEnzymeIrrevs
%modelOut.ub(ubRxns) = -model.lb(ubRxns); % in removeEnzymeIrrevs
revToIrrev = union(lbRxns, ubRxns);
modelOut.rev(revToIrrev) = 0;

% say how many bounds have changed.
REVchanged = modelOut.rev ~= model.rev;
if isfield(model, 'description')
    disp([model.description ' removed Irrevs:']);
end
disp(sprintf('Directions change: %d', length(revToIrrev))); 
