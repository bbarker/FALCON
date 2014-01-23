function generateLMoMAfluxes(model, fLabel)
% Saves a matrix of reversible-model fluxes 
% Can be used later in determining reaction direction.

nrxns = length(model.rxns);
lmomaFlux = zeros(nrxns, 60);

[celllinesarray jainMetsArray coretable] = readJainTable();
jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model);

parfor i = 1:60
    %This is copied from runComparisonScript.m; try to keep consistent:
    rxnValues = [];
    rxnList = [];
    for j = 1:length(jainMetsArray)
        jthExcIdxs = jainMetsToExcIdxs(jainMetsArray{j});
        for k = 1:length(jthExcIdxs)
            rxnList(end + 1) = jthExcIdxs(k);
            % Assume the maximum value (usually glucose) is always attainable
            % at unity:
            rxnValues(end + 1) = coretable(j, i) / max(coretable(:, i));
        end
    end
    [v_solirrev v_solrev] = runLinearMOMAOneShot(model, rxnValues, rxnList);
    lmomaFlux(:, i) = v_solrev;
end

save(['LMoMA_CoRe_Flux_' fLabel '.mat'], 'lmomaFlux')