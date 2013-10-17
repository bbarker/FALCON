function R = topCorr(Y, topN, ctype)
%
%
%INPUT
% topN   the top N variables to correlate across, given
%        that their rank is determined by the average
%        scaled absolute value across all observations.
%

% Observations correspond to different columns in Y.
%
%
if ~exist('ctype', 'var')
    ctype = 'Pearson';
end

[nr, nc] = size(Y);
% Scale the values
colSums = sum(abs(Y));
Y = median(colSums) * Y ./ repmat(colSums, nr, 1);

varMean = mean(abs(Y'));
[vMsort, vmIdx] = sort(varMean, 2, 'descend');
vMTopIdxs = vmIdx(1:topN);

R = corr(Y(vMTopIdxs, :), 'type', ctype);

