function [v_md, v_lee, v_md_s, v_lee_s, nnanDiffEnz, nnanDiffFlux, ...
         r_md, r_lee] = compareFluxByEnzymeComplexation(model,    ...
         expFile, nReps, fLabel)

%Simple script to compare fluxes output from the FALCON algorithm
%using the directly evaluation method found in Lee and the 
%min-disjunction algorithm.
%
%OUTPUT
%
%nnanDiffEnz         Number of differences in enzyme expression for
%                    the unpermuted expression data.
%
%nnanDiffFlux        Number of differences in FALCON flux for
%                    the unpermuted expression data.
%
%r_lee            Enzyme abundance esitmated by direct evaluation
%
%r_lee            Enzyme abundance estimated by FALCON.
%
%v_lee             Flux esitmated by direct evaluation.
%
%v_md             Flux estimated by FALCON.
%
%These are all saved to a .mat file (see save() below).

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

nRxns = length(model.rxns);
nRxnsIrr = length(modelIrrev.rxns);


[r_md, rs_md, r_group] = computeMinDisj(model, expFile);
[r_lee, rs_lee, r_miss] = geneToRxn(model, expFile);

[r_mdIrr, rs_mdIrr, r_groupIrr] = computeMinDisj(modelIrrev, expFile);
[r_leeIrr, rs_leeIrr, r_missIrr] = geneToRxn(modelIrrev, expFile);
%We should use the same r_group in both cases, as this is technically
%part of FALCON rather than mindisj.

[v_mdIrr, ~, ~, ~, ~, ~, v_mdIrr_s] =  ...
    falconMulti(modelIrrev, nReps, r_mdIrr, rs_mdIrr, r_groupIrr);

[v_leeIrr, ~, ~, ~, ~, ~, v_leeIrr_s] =  ...
    falconMulti(modelIrrev, nReps, r_leeIrr, rs_leeIrr, r_groupIrr);

v_md = convertIrrevFluxDistribution(v_mdIrr, matchRev);
v_lee = convertIrrevFluxDistribution(v_leeIrr, matchRev);
v_md_s = abs(convertIrrevFluxDistribution(v_mdIrr_s, matchRev));
v_lee_s = abs(convertIrrevFluxDistribution(v_leeIrr_s, matchRev));

%For testing:
%[r_md, rs_md] = deal(r_lee, rs_lee);

%Should be the same for both methods
rNotNan = ~isnan(r_md);
nnanTotal = sum(rNotNan);

r_md(isnan(r_md)) = -1;
r_lee(isnan(r_lee)) = -1;
nnanDiffEnz = sum(isDiffExp(r_md, r_lee, r_miss));
nnanDiffFlux = sum(isDiffExp(v_md, v_lee, false(nRxns, 1))); 

save(['FluxByECcomp_' model.description expFile fLabel '_' num2str(nReps) ...
      '.mat'], 'v_lee', 'v_md', 'v_md_s', 'v_lee_s', 'nnanDiffEnz',       ...
      'nnanDiffFlux', 'r_lee', 'r_lee');

% uncomment for debugging purposes:
% printDiffs(model, r_md, r_lee, r_miss);
%
% diffRxns = find(isDiffExp(r_md, r_lee, r_miss));
%dbgCell = cell(nnanDiffOrig, 4);
%for i = 1 : nnanDiffOrig
%    dbgCell{i, 1} = model.rxnNames{diffRxns(i)};
%    dbgCell{i, 2} = model.grRules{diffRxns(i)};
%    dbgCell{i, 3} = num2str(r_md(diffRxns(i)));
%    dbgCell{i, 4} = num2str(r_lee(diffRxns(i)));
%end
%cell2csv(['compareEnzymeExpression_' strrep(model.description, ' ', '_') '.csv'], dbgCell, ',');

%get expression file name

function d = isDiffExp(r1, r2, rmg)
d = columnVector(boolean(abs(r1 - r2) > 1e-4) & (~rmg))';
% end of isDiffExp

function printDiffs(m, r1, r2, rmg, priorDiffs)
% print the gene rule and computed expression levels
% (and possibly the individual gene expression levels)
% for all differences
rDiffs = find(isDiffExp(r1, r2, rmg));
if exist('priorDiffs', 'var')
    rDiffs = setdiff(rDiffs, priorDiffs);
end

for i = 1:length(rDiffs)
   ri = rDiffs(i);
   rd = num2str(abs(r1(ri) - r2(ri)));
   disp([num2str(ri) ':' 9 m.grRules{ri} 9 num2str(r1(ri)) 9 num2str(r2(ri)) 9 rd]);
end
% end of printDiffs
