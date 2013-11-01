function corrMat = compareExpressionToFlux(model, expFile)

% We need to compare across either (enzymatic) reactions or 
% genes. Since we have two datatypes related to reactions
% (estimated complex levels and fluxes), we will use the
% mean and sum of gene expression as proxies for unbiased
% gene expression levels. As the expression files should be
% preprocessed, a 0 will be assumed to be zero expression,
% not 'absent' (nan).

% since a flux can be negative, we compare the absolute
% value of the flux to the expression level.

% Enzyme complex level.


EXPCON = 0
FDEBUG = 0;
[rxn_exp_rev, rxn_exp_sd_rev, rxn_rule_group_rev] = ...
   computeMinDisj(model, expFile);

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

[rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
    computeMinDisj(modelIrrev, expFile);

v_solirrev = falcon(modelIrrev, rxn_exp, rxn_exp_sd, ...
                    rxn_rule_group, 0.01, 0, EXPCON, FDEBUG);
v_solrev = convertIrrevFluxDistribution(model, v_solirrev, matchRev);

corrType = 'Pearson';
mycorr = @(X)(corr(X, 'type', corrType));

X = [abs(v_solrev) rxn_exp_rev];

nnanRows = find(~isnan(sum(X')));
X = X(nnanRows, :);

corrMat = mycorr(X);