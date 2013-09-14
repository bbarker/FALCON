function modelConstrainedCOREExc = constrainCOREExc(model, coreValues)
%This function takes an original model, usually rec2, and outputs a cell array of models
%with the bounds on certain reactions sequentially constrained, depending on .

[celllinesarray jainMetsArray coretable] = readJainTable();
jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model);
mediumExcIdxs = loadMediumExcIdxs(model);

modelConstrainedCOREExc = model;
for i=1:length(modelConstrainedCOREExc.rxns)
    if (~isempty(regexp(modelConstrainedCOREExc.rxns{i},'^EX(.)*\(e\)$')))
        if (sum(mediumExcIdxs == i) == 0)
            modelConstrainedCOREExc.lb(i) = 0;
        end
    end
end