function [rxn_exp, rxn_exp_sd, rxn_rule_group] = ... 
    computeSimpleECexpression.m(model, genedata_filename, expFun, expStdFun)

% This function is meant to parallel computeMinDisj, but evaluates
% enzyme complex expression without taking into account gene rules.

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

