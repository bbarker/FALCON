% unfinished, but probably not needed
% would just be used to look at correlations of random permuted vectors
% given various assumptions.
function [vlens, vcorrs_mean, vcorrs_std] = randCorrTest(maxLen, step, nReps)

% Sample from a power-law distribution.
maxSZ = 1000;

vlens  = [2:step:maxLen];
vcorrs_mean = zeros(size(vlens));
vcorrs_std  = zeros(size(vlens));
parfor i = 1:length(vlens)
    corrvec = zeros(1, nReps);
    randvec = maxSZ*rand(1, vlens(i));
    for j = 1:nReps
        permVec = randperm(ngenes);
    
    end
end
