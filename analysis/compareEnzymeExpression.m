function [nnanDiffTotal nnanDiffTotalNoMiss] = ...
    compareEnzymeExpression(model, expFile, nReps)

%Simple script to compare expression similarity from
%the Lee method and the minDisj method.
%
%INPUT
%
% expFile        If equal to '', generates a random expression file.
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

deleteExpFile = false;

nRxns = length(model.rxns);
nModGenes = length(model.genes);

diffMat = zeros(nReps + 1, nRxns);
diffMatNoMiss = zeros(nReps + 1, nRxns);

if strcmp(expFile, '')
    deleteExpFile = true;
    pd = ProbDistUnivParam('gamma', [1 2]);
    gene_exp = pd.random(nModGenes, 1);
    gene_exp_sd = gene_exp .* rand(nModGenes, 1) / 10; %don't need; quick hack
    expFile = ['compareEnzymeExpressionRandData_' num2str(nReps) '.csv'];
    cellOut = cell(nModGenes + 1, 3);
    cellOut(1, :) = {'gene', 'mean', 'std'};
    cellOut(2:end, 1) = model.genes;
    cellOut(2:end, 2) = arrayfun(@num2str, gene_exp, 'UniformOutput', false);
    cellOut(2:end, 3) = arrayfun(@num2str, gene_exp_sd, 'UniformOutput', false);
    cell2csv(expFile, cellOut, '\t');
end

[r_lee, rs_lee, r_miss] = geneToRxn(model, expFile);
[r_md, rs_md, ~] = computeMinDisj(model, expFile);

%For testing:
%[r_md, rs_md] = deal(r_lee, rs_lee);

%Should be the same for both methods
rNotNan = ~isnan(r_md);
nnanTotal = sum(rNotNan);
enzTotal = sum(sum(model.rxnGeneMat') > 0);

r_md(isnan(r_md)) = -1;
r_lee(isnan(r_lee)) = -1;

diffMat(1, :) = isDiffExp(r_md, r_lee, r_miss);
diffMatNoMiss(1, :) = isDiffExpNoMiss(r_md, r_lee);

nnanDiffOrig = sum(diffMat(1, :));
nnanDiffOrigNoMiss = sum(diffMatNoMiss(1, :));

% uncomment for debugging purposes:
% printDiffs(model, r_md, r_lee, r_miss);

% diffRxns = find(isDiffExp(r_md, r_lee, r_miss));
% dbgCell = cell(nnanDiffOrig, 4);
% for i = 1 : nnanDiffOrig
%     dbgCell{i, 1} = model.rxnNames{diffRxns(i)};
%     dbgCell{i, 2} = model.grRules{diffRxns(i)};
%     dbgCell{i, 3} = num2str(r_md(diffRxns(i)));
%     dbgCell{i, 4} = num2str(r_lee(diffRxns(i)));
% end
% cell2csv(['compareEnzymeExpression_' strrep(model.description, ' ', '_') '.csv'], dbgCell, ',');

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
nnanDiffTotalNoMiss = nan*ones(1, nReps + 1);

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
        cell2csv(tmpFileName, gdout.textdata, '\t');
    elseif ndcols == 3
        %gdout.data(:,1) = genenames_perm;
        gdout.data(:,2) = gene_exp_perm;
        gdout.data(:,3) = gene_exp_sd_perm; 
        cell2csv(tmpFileName, gdout.textdata(1,:), '\t');
        dlmwrite(tmpFileName, gdout.data, '-append', 'delimiter', ...
            '\t', 'precision', 15);
    else
        disp('Problem reading gene expression file.');
        %return; ?? what to do in parfor?
    end    

    [r_lee, rs_lee, r_miss] = geneToRxn(model, tmpFileName);
    [r_md, rs_md, ~] = computeMinDisj(model, tmpFileName);

    %For testing:
    %[r_md, rs_md] = deal(r_lee, rs_lee);
    delete(tmpFileName);

    r_md(isnan(r_md)) = -1;
    r_lee(isnan(r_lee)) = -1;  
  
    % uncomment for debugging purposes:
    % (needs for not parfor loop)
    % priorDiffs = find(boolean(sum(diffMat(1 : (i-1), :))));
    % printDiffs(model, r_md, r_lee, r_miss, priorDiffs);

    diffMat(i, :) = isDiffExp(r_md, r_lee, r_miss);
    diffMatNoMiss(i, :) = isDiffExpNoMiss(r_md, r_lee);
end

nnanDiffTotal(1) = nnanDiffOrig;
nnanDiffTotalNoMiss(1) = nnanDiffOrigNoMiss;
parfor i = 2 : nReps + 1
    nnanDiffTotal(i) = sum(boolean(sum(diffMat(1 : i, :))));
    nnanDiffTotalNoMiss(i) = sum(boolean(sum(diffMatNoMiss(1 : i, :))));
end

nnanDiffAvg = mean(sum(diffMat'));
nnanDiffAvgNoMiss = mean(sum(diffMatNoMiss'));

modName = strrep(strrep(strrep(strrep(num2str(model.description), ' ', ''), ...
          '.', ''), 'xml', ''), '_', '');
[pathstr, expName, ext] = fileparts(expFile);
fileNameOut = ['compareEnzymeExpression_' modName '_' ...
                num2str(nReps) '_' expName];

disp('Saving to:');
disp(fileNameOut);
save( fileNameOut, 'nnanDiffTotal', 'nnanDiffAvg', 'nnanDiffOrig', ...
'nnanDiffTotalNoMiss', 'nnanDiffAvgNoMiss', 'nnanDiffOrigNoMiss',  ...
'nnanTotal', 'enzTotal');

if deleteExpFile %used for random files
    1; %delete(expFile); % !!!!!!!!!!!!!!!!!!!!!!!!! %
end

% end of [compareEnzymeExpression]


%apply rNotNan to diffMat

function d = isDiffExp(r1, r2, rmg)
d = columnVector(boolean(abs(r1 - r2) > 1e-4) & (~rmg))';
% end of isDiffExp

function d = isDiffExpNoMiss(r1, r2)
d = columnVector(boolean(abs(r1 - r2) > 1e-4))';
% end of isDiffExp

function printDiffs(m, r1, r2, rmg, priorDiffs)
% print the gene rule and computed expression levels
% (and possibly the individual gene expression levels)
% for all differences
rDiffs = find(isDiffExp(r1, r2, rmg));
if exist('priorDiffs', 'var')
    rDiffs = setdiff(rDiffs, priorDiffs);
end

for i = 1:length(rDiffs)
   ri = rDiffs(i);
   rd = num2str(abs(r1(ri) - r2(ri)));
   disp([num2str(ri) ':' 9 m.grRules{ri} 9 num2str(r1(ri)) 9 num2str(r2(ri)) 9 rd]);
end
% end of printDiffs
