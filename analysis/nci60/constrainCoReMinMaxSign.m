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
modelConstrained.description = [modelConstrained.description '_CoreSign'];
for k = 1:length(jainMetsArray)
    kthExcIdxs = jainMetsToExcIdxs(jainMetsArray{k});
    for j = 1:length(kthExcIdxs)
	rxnIdx = kthExcIdxs(j);
	if coretableMin(k) > 0
	    %disp(['>0: ' model.rxnNames{rxnIdx}]);
	    modelConstrained.lb(rxnIdx) = 0;
	    modelConstrained.rev(rxnIdx) = 0;
	elseif coretableMax(k) < 0
	    %disp(['<0: ' model.rxnNames{rxnIdx}]);
	    modelConstrained.ub(rxnIdx) = 0;
	    modelConstrained.rev(rxnIdx) = 0;
	end
    end
end
