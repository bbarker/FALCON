function pCorrExp = randomSeedCorrelations(model, nReps, modID, exp)
%Assumes falcon.m is using params.seed = randi...
pCorrExp = zeros(1, nReps);

parfor i = 1:nReps
    pCorrExp(i) = yeastResults(model, {'FALCON'}, exp)
end

save([modID num2str(exp) '.mat'], 'pCorrExp');