function [rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToRxn(model,genedata_filename);

% load transcript data
genedata   = importdata(genedata_filename);
gtextncols = size(genedata.textdata, 2);
gdatancols = size(genedata.data, 2);
expidx = 1;
stdidx = 2;
genenames = genedata.textdata(:, 1);
if gdatancols == 3
    genenames = ...
        arrayfun(@num2str,  genedata.data(:, 1), 'UniformOutput', false);
    expidx = 2;
    stdidx = 3;
else
    genenames(1)= [];
end

gene_exp	= genedata.data(:, expidx);
gene_exp_sd	= genedata.data(:, stdidx);

% map gene weighting to reaction weighting
[rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);
% sds 0 -> small
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;

