function modelConstrained = constrainCoReMinMaxSign(model)
% Here we check to see which exchanges reactions always have a definite
% sign across the NCI60 panel, and constrain the sign of exchange reactions
% based on this.

% Brandon Barker    09/20/2013

[celllinesarray jainMetsArray coretable] = readJainTable(1);
coretableMin = min(coretable');
coretableMax = max(coretable');
jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model);
% mediumExcIdxs = loadMediumExcIdxs(model);

modelConstrained = model;
for i = 1:length(modelConstrained.rxns)
    %If it is denoted an exchange reaction ...
    if (~isempty(regexp(modelConstrained.rxns{i}, '^EX(.)*\(e\)$')))
        foundInJainMetsToExcIdxs = 0;
	for k = 1:length(jainMetsArray)
	    kthExcIdxs = jainMetsToExcIdxs(jainMetsArray{k});
            for j = 1:length(kthExcIdxs)
                rxnIdx = kthExcIdxs(j);
                rxn = modelConstrained.rxns(rxnIdx);
                jmet = jainMetsArray{k};
		if coretableMin(k) > 0
                    modelConstrained.lb(rxnIdx) = 0;
                    % Should we unconstrain with such large bounds?
                    modelConstrained.ub(rxnIdx) = 1000;
                    modelConstrained.rev(rxnIdx) = 0;
                elseif coretableMax(k) < 0
                    modelConstrained.lb(rxnIdx) = -1000;
                    modelConstrained.ub(rxnIdx) = 0;
                    modelConstrained.rev(rxnIdx) = 0;
                end
            end
	    if(sum(kthExcIdxs == j) ~= 0)
		foundInJainMetsToExcIdxs = 1;
	    end
	end
        %For now don't constrain other reactions.
        %if (sum(mediumExcIdxs == i) == 0 && ~foundInJainMetsToExcIdxs)
        %    modelConstrained.lb(i) = 0;
        %end
    end
end