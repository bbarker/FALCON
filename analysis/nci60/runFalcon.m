function [v_solirrev v_solrev corrval] = ...
    runFalcon(model, expressionFile, rc, EXPCON, FDEBUG)
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
% rc                regularization constant on fluxes
%
% EXPCON            Flag to specify whether to use expression constraints.
%
%OUTPUT
% v_solirrev    irreversible flux vector
%
% v_solrev      reversible flux vector
%
%
% Yiping Wang    09/08/13
% Brandon Barker 09/15/2013  Now calls convertIrrevFluxDistribution

if ~exist('FDEBUG', 'var')
    FDEBUG = false;
end

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
[rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
    computeMinDisj(modelIrrev, expressionFile);

[v_solirrev corral]= falcon(modelIrrev, rxn_exp, rxn_exp_sd, ...
                    rxn_rule_group, 'rc', rc, 'minFit', 0,   ...
                    'EXPCON', EXPCON, 'FDEBUG', FDEBUG);
v_solrev = convertIrrevFluxDistribution(v_solirrev, matchRev);
