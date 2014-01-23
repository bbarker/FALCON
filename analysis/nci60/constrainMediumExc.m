function modelConstrained = constrainMediumExc(model)
% This function takes an original model, usually recon 2, and 
% constrains all models to only have very small uptake values unless
% they are one of the key medium components.

% Yiping Wang      08/??/2013
% Brandon Barker   09/21/2013    Changed to allow small amounts of
%                                uptake not listed for the medium.

minUptake = -0.005;
% The above value (-0.005) reduces growth rate, but it is still far more
% than is likely in any cancer cell line.

mediumExcIdxs = loadMediumExcIdxs(model);
modelConstrained = model;
modelConstrained.description = [modelConstrained.description '_med'];


for i = 1:length(modelConstrained.rxns)
    if (~isempty(regexp(modelConstrained.rxns{i}, '^EX(.)*\(e\)$')))
        if (sum(mediumExcIdxs == i) == 0)
            if modelConstrained.lb(i) < minUptake
                modelConstrained.lb(i) = minUptake;
            end
        end
    end
end

