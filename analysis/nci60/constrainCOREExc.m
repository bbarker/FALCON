function modelConstrained = constrainCOREExc(model)
%This function takes an original model, usually rec2, and outputs a cell array of models
%with the bounds on certain reactions sequentially constrained, depending on .

[celllinesarray jainMetsArray coretable] = readJainTable();
jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model);
mediumExcIdxs = loadMediumExcIdxs(model);

modelConstrained = model;
for i = 1:length(modelConstrained.rxns)
    if (~isempty(regexp(modelConstrained.rxns{i}, '^EX(.)*\(e\)$')))
        foundInJainMetsToExcIdxs = 0;
	for k = 1:length(jainMetsArray)
	    kthExcIdxs = jainMetsToExcIdxs(jainMetsArray{k});
	    if(sum(kthExcIdxs == i) ~= 0)
		foundInJainMetsToExcIdxs = 1;    
	    end
	end
        if (sum(mediumExcIdxs == i) == 0 && ~foundInJainMetsToExcIdxs)
            modelConstrained.lb(i) = 0;
        end
    end
end