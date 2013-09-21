function modelConstrained = testGrowthRequirements(model)
% This function is based on constrainMediumExc and checks to see
% which exchange reactions are necessary for growth (can be applied
% iteratively to get the necessary set).

% Brandon Barker    09/21/2013
mediumExcIdxs = loadMediumExcIdxs(model);
modelConstrained = model;

fbasol = optimizeCbModel(model, 'max', 'one');
for i = 1:length(modelConstrained.rxns)
    if (~isempty(regexp(modelConstrained.rxns{i}, '^EX(.)*\(e\)$')))
        if (sum(mediumExcIdxs == i) == 0)
            modelConstrained.lb(i) = 0;
            if fbasol.x(i) < 0
                disp(modelConstrained.rxnNames{i});
                fbasol = optimizeCbModel(modelConstrained, 'max', 'one');
                disp(['Reaction: ' num2str(i)]);
                disp([modelConstrained.lb(i), ... 
                      modelConstrained.ub(i) fbasol.f]);
                if fbasol.f < 1e-7
                   break
                end
            end
    end
end