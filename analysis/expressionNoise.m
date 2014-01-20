function [sigmaVec, PcorrV, ScorrV, KcorrV, PcorrE, ScorrE, KcorrE, ...
          PcorrEV, ScorrEV, KcorrEV, Vactive, VdeltaSign, VOnOff,   ...
          VnormDiff, EnormDiff, TimeRec, IterRec] =                          ...
          expressionNoise(model, expFile, sigMax, reps, simLabel, LPmeth)

regC   = 0;
expCon = false;

if ~exist('simLabel', 'var')
    simLabel = '';
end
if ~exist('LPmeth', 'var')
    LPmeth = 1;
end

LPseed = 0; % Use default seed; don't need stochasticity in solver since we
            % are testing stochasticity due to expression noise.

% Things to track:
%1: Correlation between e and e_perturbed
%2: Noise level e.g. sigma for gaussian error
%  For now do sigma = [0.01 0.05 0.1 0.5 1 1.5 2] ?
%  (Consider covariance effects later?  For now assume independent multiplicative error on expression) 
%3: Number of "non-zero" fluxes.
%4: Number of flux changes (positive <-> negative and non-zero <-> zero)
%5: Correlation between v and v_perturbed (pearson, spearman, kendall)
%6: norm diff between e and e_p, v and v_p

% TODO
% Randomly select sigma to add more variation
% Confidence intervals? see jason's email
% Fix to work with Lung data

sigmaVec = sigMax * rand(1, reps);

nSigmas = length(sigmaVec);
PcorrV = zeros(1, reps);
ScorrV = zeros(1, reps);
KcorrV = zeros(1, reps);
PcorrE = zeros(1, reps);
ScorrE = zeros(1, reps);
KcorrE = zeros(1, reps);
PcorrEV = zeros(1, reps);
ScorrEV = zeros(1, reps);
KcorrEV = zeros(1, reps);
Vactive = zeros(1, reps);
VdeltaSign = zeros(1, reps);
VOnOff  = zeros(1, reps);
VnormDiff  = zeros(1, reps);
EnormDiff  = zeros(1, reps);
TimeRec    = zeros(1, reps);
IterRec    = zeros(1, reps);

%Make irreversible model
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
nIrxns = length(modelIrrev.rxns);
[rxn_exp, rxn_exp_sd, rxn_rule_group] = computeMinDisj(modelIrrev, expFile);

%errMat = zeros(nSigmas, reps, nIrxns);
%parfor i = 1:nSigmas
%  errMat(i,:,:) = reshape(simpleTruncatedNorm(sigmaVec(i), 0, inf, reps*nIrxns, 1), reps, nIrxns);
%end



%Now compute initial flux:
[v, corrval] = falcon(modelIrrev, rxn_exp, rxn_exp_sd, rxn_rule_group, ...
                      'rc', regC, 'EXPCON', expCon,  ...
                      'LPmeth', LPmeth, 'LPseed', LPseed); 
if norm(v,1) < 1e-7
    disp('Error, initial flux prediction failed.');
    return;
end

vrev = convertIrrevFluxDistribution(v, matchRev);

% = simpleTruncatedNorm(sigma, a, b, , mu)
Zthresh = 1e-6;

vrevSign = signThresh(vrev, Zthresh);

parfor i = 1:reps
    [rxn_exp_p, ~, ~] = computeMinDisj(modelIrrev, expFile,  sigmaVec(i));
    [rxn_exp_rev_p, ~, ~] = computeMinDisj(model, expFile,  sigmaVec(i));

    %It seems there is a small difference in the number of pertubed nans from time to time, so
    %we need to check for the common not-nans. However, it is not clear to me why the perturbation
    %can cause this, so it should be checked further.
    e_not_nan = boolean((~isnan(rxn_exp_p)) .* (~isnan(rxn_exp)));
    PcorrE(i) = corr(rxn_exp(e_not_nan), rxn_exp_p(e_not_nan), 'type', 'Pearson');
    ScorrE(i) = corr(rxn_exp(e_not_nan), rxn_exp_p(e_not_nan), 'type', 'Spearman'); 
    KcorrE(i) = corr(rxn_exp(e_not_nan), rxn_exp_p(e_not_nan), 'type', 'Kendall');
    EnormDiff(i) = norm(rxn_exp(~isnan(rxn_exp + rxn_exp_p)) - ...
                   rxn_exp_p(~isnan(rxn_exp_p + rxn_exp)), 1);

    [v_p, corrval_p, n, VA, fTime, fIter] = falcon(modelIrrev, rxn_exp_p, ... 
        rxn_exp_sd, rxn_rule_group, 'rc', regC, 'EXPCON', expCon, ...
        'LPmeth', LPmeth, 'LPseed', LPseed);
    TimeRec(i) = fTime;
    IterRec(i) = fIter;
    vrev_p = convertIrrevFluxDistribution(v_p, matchRev);
    PcorrV(i) = corr(vrev(:), vrev_p(:), 'type', 'Pearson');
    ScorrV(i) = corr(vrev(:), vrev_p(:), 'type', 'Spearman'); 
    KcorrV(i) = corr(vrev(:), vrev_p(:), 'type', 'Kendall');

    PcorrEV(i) = corr(rxn_exp_rev_p(~isnan(rxn_exp_rev_p)), ...
        abs(vrev_p(~isnan(rxn_exp_rev_p))), 'type', 'Pearson');
    ScorrEV(i) = corr(rxn_exp_rev_p(~isnan(rxn_exp_rev_p)), ...
        abs(vrev_p(~isnan(rxn_exp_rev_p))), 'type', 'Spearman'); 
    KcorrEV(i) = corr(rxn_exp_rev_p(~isnan(rxn_exp_rev_p)), ...
        abs(vrev_p(~isnan(rxn_exp_rev_p))), 'type', 'Kendall');
    Vactive(i) = sum(abs(vrev_p)>Zthresh);
    vrevSign_p = signThresh(vrev_p, Zthresh);
    VdeltaSign(i) = sum(vrevSign_p .* vrevSign < 0);
    VOnOff(i) = sum((abs(vrevSign_p) + abs(vrevSign)) == 1);
end

finishTime = strrep(strrep(num2str(clock()), ' ', ''), '.', '');
modName = strrep(strrep(strrep(strrep(num2str(model.description), ' ', ''), ...
          '.', ''), 'xml', ''), '_', '');
[pathstr, expName, ext] = fileparts(expFile);

save(['pertData_' modName '_' expName '_' num2str(sigMax) '_' num2str(reps) ...
      '_' finishTime '_' simLabel '_' num2str(LPmeth) '.mat'], 'sigmaVec',  ...
     'PcorrV', 'ScorrV', 'KcorrV', 'PcorrE', 'ScorrE', 'KcorrE', 'PcorrEV', ... 
     'ScorrEV', 'KcorrEV', 'Vactive', 'VdeltaSign', 'VOnOff', 'VnormDiff',  ...
     'EnormDiff', 'TimeRec', 'IterRec');  

