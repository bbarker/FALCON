function [v_sol, lIter, v_sol_s, lIter_s, lTime] = ...
    dataToFluxFixMulti(m, r, r_sd, nReps)

% This is like dataToFlux, but returns a mean and std
% for v_sol and lIter, and saves the flux distributions.

nRxns = length(m.rxns);


v_sol_Dist   = zeros(nReps, nRxns);
lIter_Dist   = zeros(nReps, 1);
lTimeVec     = zeros(nReps, 1);

timeInit = num2str(now());

%normally parfor, disabled for testing
for i = 1:nReps
    disp(['Starting dataToFlux rep: ' num2str(i)]);
    dFstart = tic;
    [v_sol, lIter]    = dataToFluxFix(m, r, r_sd);
    lTimeVec(i) = toc(dFstart);
    v_sol_Dist(i, :)  = columnVector(v_sol)';
    lIter_Dist(i)     = lIter;
end

v_sol   = mean(v_sol_Dist)';
lIter   = mean(lIter_Dist);

v_sol_s   = std(v_sol_Dist)';
lIter_s   = std(lIter_Dist);
lTime = mean(lTimeVec);

fileNameOut = ['dataToFlux_' num2str(nReps) '_' timeInit '.mat'];
save(fileNameOut, 'v_sol', 'v_sol_s', 'v_sol_Dist', ...
    'lIter', 'lIter_s', 'lIter_Dist', 'lTimeVec', 'lTime');

