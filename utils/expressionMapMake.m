function [gMap, sMap] = expressionMapMake(model, expFile)
% WARNING: untested with NaN values

% This function simply returns a Map container with
% keys == gene names
% values == expression levels

genedata	= importdata(expFile);
genenames	= genedata.textdata(:, 1);
genenames(1)= [];
gene_exp	= genedata.data(:,1);
gene_exp_sd	= genedata.data(:,2);

gMap = containers.Map();
sMap = containers.Map();
for i = 1:length(genenames)
    gMap(genenames{i}) = gene_exp(i);
    sMap(genenames{i}) = gene_exp_sd(i);
end
