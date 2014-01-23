function modelConstrained = constrainImputedInternal(model, allRxns, LMoMA_file, expFile)
% These constraints are used to constrain the direction of internal,
% reactions based on some form of imputation (for now, just linear MoMA).

% We constrain only reaction directions consistent across all cell lines,

%OPTIONAL INPUTS
% allRxns    If allRxns is set, we constrain not only the sign of  
%            enzymatic reactions, but any reaction with uniform flux
% expFile    Expression file that is used with computeMinDisj to 
%            see how many constrained reactions actually have
%            available data in a given experimental dataset.

if ~exist('allRxns', 'var')
    allRxns = false;
end

%Load LMoMA flux (see generateLMoMAfluxes.m)
%load('LMoMA_CoRe_Flux.mat', 'lmomaFlux');
load(LMoMA_file, 'lmomaFlux');

nrxns = length(model.rxns);
[selExc, selUpt] = findExcRxns(model);
geneRxns = boolean(sum(model.rxnGeneMat')');
consSign = sign(lmomaFlux(:, 1)) .* (all(lmomaFlux' < 0) | all(lmomaFlux' > 0))';
alllb = '';
if ~allRxns
    consSign = consSign .* ~selExc .* geneRxns;
    alllb = '_not_all';
else
    alllb = '_all';
end

nNegDir = full(sum(consSign < 0))
nPosDir = full(sum(consSign > 0))

modelConstrained = model;
modelConstrained.description = [modelConstrained.description '_lmoma' alllb];
for i = 1:nrxns
    if consSign(i) < 0
        modelConstrained.ub(i) = 0;
        modelConstrained.rev(i) = false;
    elseif consSign(i) > 0
        modelConstrained.lb(i) = 0;
        modelConstrained.rev(i) = false;
    end
end

if exist('expFile', 'var')
    % I'll take it as a sign you want a report.
    [pathstr, expName, ext] = fileparts(expFile);
    [rxn_exp, rxn_exp_sd, rxn_rule_group] = computeMinDisj(model, expFile);

    arLabel = '';
    if allRxns
        arLabel = 'allReactions';
    end
    reportFI = fopen(['LMoMA_CoRe_report' arLabel '_' expName '.csv'], 'w');
    fprintf(reportFI, '%s\t%d\n', '#Forward reactions inferred:', nPosDir);
    fprintf(reportFI, '%s\t%d\n', '#Backward reactions inferred:', nNegDir);
    
    consSignExp = consSign .* sign(rxn_exp); 
    nNegDirExp = full(sum(consSignExp < 0))
    nPosDirExp = full(sum(consSignExp > 0))
    fprintf(reportFI, '%s\t%d\n', ...
        '#Forward reactions inferred with expression:', nPosDirExp);
    fprintf(reportFI, '%s\t%d\n', ...
        '#Backward reactions inferred with expression:', nNegDirExp);

    fprintf(reportFI, '%s\t%d\n', '#Bounds changed in model:', ...
            sum((model.lb ~= modelConstrained.lb) | ...
                (model.ub ~= modelConstrained.ub)));
    fprintf(reportFI, '\n%s\t%s\t%s\n', 'reaction name', ...
        'direction', 'EC intensity');
    for i = 1:nrxns
        if consSign(i)
            fprintf(reportFI, '%s\t%f\t%f\n', model.rxnNames{i}, ...
                full(consSign(i)), rxn_exp(i));
        end
    end
end

