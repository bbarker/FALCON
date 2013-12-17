function [nnanDiffTotal, nnanDiffAvg, nnanDiffOrig, nnanTotal] = ...
    compareEnzymeExpression(model, expFile, nReps)
testx = 1;
%Simple script to compare expression similarity from
%the Lee method and the minDisj method.
%
%OUTPUT
%nnanDiffTotal    Total number of differences in expression
%                 found across all permutations of enzyme
%                 expression, not counting multiple differences
%                 occurring for the same reaction.
%
%nnanDiffAvg      Average number of differences in expression
%                 across all permutations of enzyme expression.
%
%nnanDiffOrig     Number of differences in enzyme expression for
%                 the unpermuted expression data.
%
%nnanTotal        For reference, the total number of reactions
%                 in the unpermuted condition with expression data.

nRxns = length(model.rxns);

diffMat = zeros(nReps + 1, nRxns);

[r_lee, rs_lee, ~] = geneToRxn(model, expFile);
[r_md, rs_md, ~] = computeMinDisj(model, expFile);


%Should be the same for both methods
rNotNan = ~isnan(r_md);
nnanTotal = sum(rNotNan);

nnanDiffOrig = sum(r_md(rNotNan) ~= r_lee(rNotNan));
diffMat(1, :) = (r_md ~= r_lee)';

for i = 2 : nReps + 1
    nnanDiffAvg = 0;
    nnanDiffTotal = 0;
end

%apply rNotNan to diffMat