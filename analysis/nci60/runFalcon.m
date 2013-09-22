function [v_solirrev v_solrev numiterations] = ...
    runFalcon(model, expressionFile, flux_sum, rc)
%INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   lb           Lower bounds
%   ub           Upper bounds
%   grRules      boolean rules in terms of genes (see next field)
%   genes        Entrez gene ids (without isoform id)
%
% expressionFile    Tab-delimited file with a head for the columns:
%                   gene (entrez gene id), mean (expression value,
%                   and standard deviation (of expression).
%
% flux_sum          Positive value used in the FALCON objective
%                   to guarantee a non-zero flux vector with 
%                   1-norm exceeding flux_sum is chosen.
%
% rc                regularization constant on fluxes
%
%OUTPUT
% v_solirrev    irreversible flux vector
%
% v_solrev      reversible flux vector
%
%
% Yiping Wang    09/08/13
% Brandon Barker 09/15/2013  Now calls convertIrrevFluxDistribution

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
[rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
    computeMinDisj(modelIrrev, expressionFile);

v_solirrev = falcon(modelIrrev, rxn_exp, rxn_exp_sd, ...
                    rxn_rule_group, flux_sum, rc, 0);
v_solrev = convertIrrevFluxDistribution(model, v_solirrev, matchRev);
