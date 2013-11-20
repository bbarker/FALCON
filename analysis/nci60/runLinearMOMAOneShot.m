function [v_solirrev v_solrev] = ...
    runLinearMOMAOneShot(model, rxnValues, rxnList)
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
WTflux = zeros(length(model.rxns),1);
rxnListIrrev = [];
for i = 1:length(rxnList)
    WTflux(rxnList(i)) = rxnValues(i);
    rxnListIrrev = [rxnListIrrev rev2irrev{rxnList(i)}];
end

WTflux = convertToIrreversibleFluxVector(modelIrrev, irrev2rev, WTflux);

[solutionDel, solutionWT] = linearMOMAOneShot(modelIrrev, WTflux, ...
                                0, 0, 1, ...
                                zeros(length(modelIrrev.rxns), 1), ...
                                rxnListIrrev);
v_solirrev = solutionDel.x;
v_solrev = convertIrrevFluxDistribution(v_solirrev, matchRev);