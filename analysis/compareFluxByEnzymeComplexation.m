function [Fdirect, Ffalcon, nnanDiffEnz, nnanDiffFlux,  ...
          ECdirect, ECfalcon] = ...
          compareFluxByEnzymeComplexation(model, expFile, fLabel)

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
%ECdirect            Enzyme abundance esitmated by direct evaluation
%
%ECfalcon            Enzyme abundance estimated by FALCON.
%
%Fdirect             Flux esitmated by direct evaluation.
%
%Ffalcon             Flux estimated by FALCON.
%
%These are all saved to a .mat file (see save() below).

nRxns = length(model.rxns);

[r_lee, rs_lee, r_miss] = geneToRxn(model, expFile);
[r_md, rs_md, ~] = computeMinDisj(model, expFile);


if useMinDisj % need to modify to use old expression values
    [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s] = ...
	falconMulti(modelIrrev, rxn_exp_md, ...
	rxn_exp_sd_md, rxn_rule_group, nReps, regC, minFit, expCon);
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);

else
    [v_directIrr, ~, ~, ~, ~, ~, v_directIrr_s] = ...
	falconMulti(modelIrrev, rxn_exp_irr, rxn_exp_sd_irr, ...
	rxn_rule_group, nReps, regC, minFit, expCon);
    v_direct = convertIrrevFluxDistribution(v_directIrr, matchRev);

end



%For testing:
%[r_md, rs_md] = deal(r_lee, rs_lee);

%Should be the same for both methods
rNotNan = ~isnan(r_md);
nnanTotal = sum(rNotNan);





r_md(isnan(r_md)) = -1;
r_lee(isnan(r_lee)) = -1;
nnanDiffEnz = sum(isDiffExp(r_md, r_lee, r_miss));

save(['FluxByECcomp_' model.description expFile fLabel '.mat'], ...
    'Fdirect', 'Ffalcon', 'nnanDiffEnz', 'nnanDiffFlux', ...
    'ECdirect', 'ECfalcon');


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
%cell2csv(['compareEnzymeExpression_' strrep(model.description, ' ', '_') '.csv'], dbgCell, ',', 2000);

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
