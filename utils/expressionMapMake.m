function [getGeneExp, getGeneVar] = expressionMapMake(expFile)

% This function simply returns a Map container with
% keys == gene names
% values == expression levels

genedata         = importdata(expFile);
[ndrows, ndcols] = size(genedata.data);
if ndcols == 2
    genenames	= genedata.textdata(:, 1);
    genenames(1)= [];
    gene_exp	= genedata.data(:,1);
    gene_exp_sd	= genedata.data(:,2);
elseif ndcols == 3
    genenames	= cellfun(@num2str, ...
        columnVector(num2cell(genedata.data(:, 1))), 'UniformOutput', false);
    gene_exp	= genedata.data(:,2);
    gene_exp_sd	= genedata.data(:,3);
end

gMap = containers.Map();
vMap = containers.Map();
for i = 1:length(genenames)
    gMap(genenames{i}) = gene_exp(i);
    vMap(genenames{i}) = gene_exp_sd(i);
end
gMap('nan') = nan;
gMap('NaN') = nan;
gMap('NAN') = nan;
vMap('nan') = nan;
vMap('NaN') = nan;
vMap('NAN') = nan;

function gV = makeGetGeneVal(aMap)
    function exp = getGeneVal(key)
        if aMap.isKey(key) 
            exp = aMap(key); 
        else 
            exp = nan;
        end
    end
    gV = @getGeneVal;
end % of makeGetGeneVal

getGeneExp = makeGetGeneVal(gMap);
getGeneVar = makeGetGeneVal(vMap);

end % of expressionMapMake



