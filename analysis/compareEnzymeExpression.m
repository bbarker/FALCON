function [nnanDiffTotal, nnanDiffAvg, nnanDiffOrig, nnanTotal] = ...
    compareEnzymeExpression(model, expFile, nReps)
testx = 1;
%Simple script to compare expression similarity from
%the Lee method and the minDisj method.
%
%OUTPUT
%nnanDiffTotal    Total number of differences in expression
%                 found across all permutations of enzyme
%                 expression, not counting multiple differences
%                 occurring for the same reaction.
%
%nnanDiffAvg      Average number of differences in expression
%                 across all permutations of enzyme expression.
%
%nnanDiffOrig     Number of differences in enzyme expression for
%                 the unpermuted expression data.
%
%nnanTotal        For reference, the total number of reactions
%                 in the unpermuted condition with expression data.

nRxns = length(model.rxns);

diffMat = zeros(nReps + 1, nRxns);

[r_lee, rs_lee, ~] = geneToRxn(model, expFile);
[r_md, rs_md, ~] = computeMinDisj(model, expFile);
%For testing:
%[r_md, rs_md] = deal(r_lee, rs_lee);

%Should be the same for both methods
rNotNan = ~isnan(r_md);
nnanTotal = sum(rNotNan);

r_md(isnan(r_md)) = -1;
r_lee(isnan(r_lee)) = -1;


nnanDiffOrig = sum(r_md ~= r_lee);
diffMat(1, :) = (r_md ~= r_lee)';

% uncomment for debugging purposes:
printDiffs(model, r_md, r_lee);

diffRxns = find(r_md ~= r_lee);
dbgCell = cell(nnanDiffOrig, 4);
for i = 1 : nnanDiffOrig
    dbgCell{i, 1} = model.rxnNames{diffRxns(i)};
    dbgCell{i, 2} = model.grRules{diffRxns(i)};
    dbgCell{i, 3} = num2str(r_md(diffRxns(i)));
    dbgCell{i, 4} = num2str(r_lee(diffRxns(i)));
end
cell2csv(['compareEnzymeExpression_' strrep(model.description, ' ', '_') '.csv'], dbgCell, ',', 2000);

%get expression file name
[pathstr, expName, ext] = fileparts(expFile);

genedata         = importdata(expFile);
[ndrows, ndcols] = size(genedata.data);
if ndcols == 2
    genenames	= genedata.textdata(:,1);
    genenames(1)= [];
    gene_exp	= genedata.data(:,1);
    gene_exp_sd	= genedata.data(:,2);
elseif ndcols == 3
    %genenames	= cellfun(@num2str, ...
    %    columnVector(num2cell(genedata.data(:, 1))), 'UniformOutput', false);
    genenames    = genedata.data(:,1); 
    %genenames(1)= [];
    gene_exp	 = genedata.data(:,2);
    gene_exp_sd	 = genedata.data(:,3);
end

ngenes = length(gene_exp);
nnanDiffTotal = nan*ones(1, nReps + 1);

parfor i = 2 : (nReps + 1)
    permVec = randperm(ngenes);
    gene_exp_perm = gene_exp(permVec);
    gene_exp_sd_perm = gene_exp_sd(permVec);
    %genenames_perm = genenames(permVec);
    % Since FALCON currently only supports calling minDisj on a file,
    % this means we need to write out all the data to a randomly
    % generated file and then delete the file.
    tmpFileName = ['compareEnzymeExpression_' num2str(nReps) '_' num2str(i) ...
                   '_' expName '.csv'];
    gdout = genedata;
    if ndcols == 2
        % need permuted gene names too
        % gdout.textdata(2:end, 1) = genenames_perm;
        gdout.textdata(2:end, 2) = cellfun(@num2str, ...
            num2cell(gene_exp_perm), 'UniformOutput', false);
        gdout.textdata(2:end, 3) = cellfun(@num2str, ...
            num2cell(gene_exp_sd_perm), 'UniformOutput', false);
        cell2csv(tmpFileName, gdout.textdata, '\t', 2000);
    elseif ndcols == 3
        %gdout.data(:,1) = genenames_perm;
        gdout.data(:,2) = gene_exp_perm;
        gdout.data(:,3) = gene_exp_sd_perm; 
        cell2csv(tmpFileName, gdout.textdata(1,:), '\t', 2000);
        dlmwrite(tmpFileName, gdout.data, '-append', 'delimiter', ...
            '\t', 'precision', 15);
    else
        disp('Problem reading gene expression file.');
        %return; ?? what to do in parfor?
    end    

    [r_lee, rs_lee, ~] = geneToRxn(model, tmpFileName);
    [r_md, rs_md, ~] = computeMinDisj(model, tmpFileName);
    %For testing:
    %[r_md, rs_md] = deal(r_lee, rs_lee);
    delete(tmpFileName);

    r_md(isnan(r_md)) = -1;
    r_lee(isnan(r_lee)) = -1;    
    diffMat(i, :) = (r_md ~= r_lee)';
end

nnanDiffTotal(1) = nnanDiffOrig;
parfor i = 2 : nReps + 1
    nnanDiffTotal(i) = sum(boolean(sum(diffMat(1 : i, :))));
end

nnanDiffAvg = mean(sum(diffMat'));

%apply rNotNan to diffMat

function printDiffs(m, r1, r2)
% print the gene rule and computed expression levels
% (and possibly the individual gene expression levels)
% for all differences

rDiffs = find(r1 ~= r_2);

for i = 1:length(rDiffs)
   ri = rDiffs(i);
   disp([m.grRules(ri) '\t' num2str(r1(ri)) '\t' num2str(r2(ri))]);
end
