function [rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToRxn(model,genedata_filename);

% load transcript data
genedata	= importdata(genedata_filename);
genenames	= genedata.textdata(:,1);
genenames(1)= [];
gene_exp	= genedata.data(:,1);
gene_exp_sd	= genedata.data(:,2);

% map gene weighting to reaction weighting
[rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);
% sds 0 -> small
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;

