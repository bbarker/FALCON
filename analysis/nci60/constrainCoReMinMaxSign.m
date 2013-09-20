function modelConstrained = constrainCoReMinMaxSign(model)
% Here we check to see which exchanges reactions always have a definite
% sign across the NCI60 panel, and constrain the sign of exchange reactions
% based on this.

% Brandon Barker    09/20/2013

[celllinesarray jainMetsArray coretable] = readJainTable();
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
	    kthExcIdxs = jainMetsToExcIdxs(jainMetsArray{k})
	    if(sum(kthExcIdxs == j) ~= 0)
		foundInJainMetsToExcIdxs = 1;
	    end
	end
%	disp (foundInJainMetsToExcIdxs)
        if foundInJainMetsToExcIdxs
	    disp(size(find(kthExcIdxs)))
	    
            %modelConstrained.lb(i)
            %modelConstrained.ub(i)
        end
        %For now don't constrain other reactions.
        %if (sum(mediumExcIdxs == i) == 0 && ~foundInJainMetsToExcIdxs)
        %    modelConstrained.lb(i) = 0;
        %end
    end
end