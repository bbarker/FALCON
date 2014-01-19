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


EXPCON = false
FDEBUG = false;
[rxn_exp_rev, rxn_exp_sd_rev, rxn_rule_group_rev] = ...
   computeMinDisj(model, expFile);

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

[rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
    computeMinDisj(modelIrrev, expFile);

v_solirrev = falcon(modelIrrev, rxn_exp, rxn_exp_sd, rxn_rule_group, ...
     'rc', 0.01, 'EXPCON', EXPCON, 'FDEBUG', FDEBUG);

v_solrev = convertIrrevFluxDistribution(v_solirrev, matchRev);

function sumout = mysum(x)
    if numel(x) == 0
        sumout = nan;
    else
        sumout = sum(x);
    end
end

function minout = mymin(x)
    if numel(x) == 0
        minout = nan;
    else
        minout = min(x);
    end
end

function maxout = mymax(x)
    if numel(x) == 0
        maxout = nan;
    else
        maxout = max(x);
    end
end

[rxn_exp_mean, rxn_exp_sd_mean] = computeSimpleECexpression(...
    model, expFile, @mean, @mean);

[rxn_exp_sum, rxn_exp_sd_sum] = computeSimpleECexpression(...
    model, expFile, @mysum, @mysum);

[rxn_exp_min, rxn_exp_sd_min] = computeSimpleECexpression(...
    model, expFile, @mymin, @mymin);

[rxn_exp_max, rxn_exp_sd_max] = computeSimpleECexpression(...
    model, expFile, @mymax, @mymax);

[rxn_exp_median, rxn_exp_sd_median] = computeSimpleECexpression(...
    model, expFile, @median, @median);

corrType = 'Pearson';
mycorr = @(X)(corr(X, 'type', corrType));

X = [abs(v_solrev) rxn_exp_rev rxn_exp_sum rxn_exp_mean rxn_exp_median ...
     rxn_exp_min rxn_exp_max];

save('Xcompare.mat', 'X');

nnanRows = find(~isnan(sum(X')));
X = X(nnanRows, :);

sz_X = size(X)

corrMat = mycorr(X);

end