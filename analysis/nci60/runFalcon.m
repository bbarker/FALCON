function [v_solirrev v_solrev numiterations]=runFalcon(model,expressionFile,flux_sum,rc)
% Yiping Wang    09/08/13
% Brandon Barker 09/15/2013  Now calls convertIrrevFluxDistribution
[modelIrrev,matchRev,rev2irrev,irrev2rev]=convertToIrreversible(model);
[rxn_exp,rxn_exp_sd,rxn_rule_group]=computeMinDisj(modelIrrev,expressionFile);

v_solirrev=falcon(modelIrrev,rxn_exp,rxn_exp_sd,rxn_rule_group,flux_sum,rc,0);
v_solrev = convertIrrevFluxDistribution(model, v_solirrev, matchRev);
