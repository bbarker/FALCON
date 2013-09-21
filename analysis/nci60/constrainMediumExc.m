function modelConstrained = constrainMediumExc(model)
%This function takes an original model, usually rec2, and outputs a cell array of models
%with the bounds on certain reactions sequentially constrained, depending on .

mediumExcIdxs = loadMediumExcIdxs(model);

modelConstrained = model;
for i = 1:length(modelConstrained.rxns)
    if (~isempty(regexp(modelConstrained.rxns{i}, '^EX(.)*\(e\)$')))
        if (sum(mediumExcIdxs == i) == 0)
            modelConstrained.lb(i) = 0;
        end
    end
end