function [v_solirrev v_solrev] = runFBA(model)
%INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%
%OUTPUT
% v_solirrev    irreversible flux vector
%
% v_solrev      reversible flux vector
%
%
% Yiping Wang    09/08/13
% Brandon Barker 09/14/2013  Now calls convertIrrevFluxDistribution

[modelIrrev matchRev rev2irrev irrev2rev] = convertToIrreversible(model);
modelIrrev = changeObjective(modelIrrev, 'biomass_reaction');
solutionStructIrrev = optimizeCbModel(modelIrrev, 'max');
v_solirrev = solutionStructIrrev.x;
v_solrev = convertIrrevFluxDistribution(v_solirrev, matchRev);