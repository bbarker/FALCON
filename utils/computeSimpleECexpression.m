function [rxn_exp, rxn_exp_sd] = ... 
    computeSimpleECexpression(model, genedata_filename, expFun, expVarFun)

% This function is meant to parallel computeMinDisj, but evaluates
% enzyme complex expression without taking into account gene rules.
% rxn_rule_group can just be obtained from running computeMinDisj.m
% 
%INPUT
% model                 A cobra model with GPRs
%
% genedata_filename     Expression file with mean and stdevs.
%
% expFun                A function to evaluate the enzyme complex expression
%                       given an array of gene expression levels
%                       associated to the complex.
%                      
% expStdFun             A function to evaluate the enzyme complex expression
%                       standard deviation given an array of gene expression 
%                       standard deviations associated to the complex.
%

%Note, this implementation is currently quite slow, but that is ok as we 
%are only intending on using this for some minimal testing purposes.

[getGeneExp, getGeneVar] = expressionMapMake(genedata_filename);

nrxns = length(model.rxns);
rxn_exp = zeros(nrxns, 1);
rxn_exp_sd = zeros(nrxns, 1);

for i = 1:nrxns
    rxnGeneIdxs = find(model.rxnGeneMat(i, :));
    rxn_exp(i) = expFun(cellfun(getGeneExp, model.genes(rxnGeneIdxs)));
    rxn_exp_sd(i) = expVarFun(cellfun(getGeneVar, model.genes(rxnGeneIdxs)));
end

