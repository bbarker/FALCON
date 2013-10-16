function [numCovered, percCovered] = checkReactionsCovered(model, expressionFile)

numCovered = 0;
percCovered = 0;

%[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
[rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
    computeMinDisj(model, expressionFile);

ruleLen = cellfun(@length, model.grRules);
nRuleRxns = sum(ruleLen > 0);
numCovered = sum(~isnan(rxn_exp));
percCovered = 100*numCovered/nRuleRxns;