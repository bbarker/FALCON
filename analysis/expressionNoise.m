function [sigmaOut, PcorrV, ScorrV, KcorrV, PcorrE, ScorrE, KcorrE, PcorrEV, ScorrEV, KcorrEV, ...
	  Vactive, VdeltaSign, VOnOff, VnormDiff, EnormDiff] = expressionNoise(model, expFile, sigmaVec, reps)
rxnScale = {'pyruvate kinase'}

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

nSigmas = length(sigmaVec);
sigmaOut = zeros(1, nSigmas*reps); 
PcorrV = zeros(1, nSigmas*reps);
ScorrV = zeros(1, nSigmas*reps);
KcorrV = zeros(1, nSigmas*reps);
PcorrE = zeros(1, nSigmas*reps);
ScorrE = zeros(1, nSigmas*reps);
KcorrE = zeros(1, nSigmas*reps);
PcorrEV = zeros(1, nSigmas*reps);
ScorrEV = zeros(1, nSigmas*reps);
KcorrEV = zeros(1, nSigmas*reps);
Vactive = zeros(1, nSigmas*reps);
VdeltaSign = zeros(1, nSigmas*reps);
VOnOff  = zeros(1, nSigmas*reps);
VnormDiff  = zeros(1, nSigmas*reps);
EnormDiff  = zeros(1, nSigmas*reps);



%Add a check (especially for initial flux) to skip if flux is zero and generate a new randomization
%Also, report an error (count)

nSims = nSigmas*reps;


%Make irreversible model
[modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);
nIrxns = length(modelIrrev.rxns);
[rxn_exp,rxn_exp_sd,rxn_rule_group] = computeMinDisj(modelIrrev, expFile);
%errMat = zeros(nSigmas, reps, nIrxns);
%parfor i = 1:nSigmas
%  errMat(i,:,:) = reshape(simpleTruncatedNorm(sigmaVec(i), 0, inf, reps*nIrxns, 1), reps, nIrxns);
%end

%Now compute initial flux:
[v, corrval, nvar] = eMOMA6(modelIrrev,rxn_exp,rxn_exp_sd,rxn_rule_group,rxnScale,0,{1});
vrev = convertIrrevFluxDistribution(v,matchRev);

% = simpleTruncatedNorm(sigma, a, b, , mu)
Zthresh = 1e-6;

vrevSign = signThresh(vrev, Zthresh);

parfor ij = 1:nSims
  ijp = ij-1;  % ijp in [0, 99]
  sigmaIdx = mod(ijp,nSigmas) + 1;
  rep = floor(ijp/nSigmas) + 1;   
  [rxn_exp_p,~,~] = computeMinDisj(modelIrrev, expFile,  sigmaVec(sigmaIdx));
  [rxn_exp_rev_p,~,~] = computeMinDisj(model, expFile,  sigmaVec(sigmaIdx));

  %It seems there is a small difference in the number of pertubed nans from time to time, so
  %we need to check for the common not-nans. However, it is not clear to me why the perturbation
  %can cause this, so it should be checked further.
  e_not_nan = boolean((~isnan(rxn_exp_p)) .* (~isnan(rxn_exp)));
  disp('rxn_exp sizes');
  disp([sigmaVec(sigmaIdx) length(rxn_exp) sum(isnan(rxn_exp)) length(rxn_exp_p) sum(isnan(rxn_exp_p))]);
  PcorrE(ij) = corr(rxn_exp(e_not_nan), rxn_exp_p(e_not_nan), 'type', 'Pearson');
  ScorrE(ij) = corr(rxn_exp(e_not_nan), rxn_exp_p(e_not_nan), 'type', 'Spearman'); 
  KcorrE(ij) = corr(rxn_exp(e_not_nan), rxn_exp_p(e_not_nan), 'type', 'Kendall');
  EnormDiff(ij) = norm(rxn_exp(~isnan(rxn_exp)) - rxn_exp_p(~isnan(rxn_exp_p)), 1);
  sigmaOut(ij) = sigmaVec(sigmaIdx);

  [v_p, corrval_p, nvar_p] = eMOMA6(modelIrrev,rxn_exp_p,rxn_exp_sd,rxn_rule_group,rxnScale,0,{1});
  vrev_p = convertIrrevFluxDistribution(v_p,matchRev);
  PcorrV(ij) = corr(vrev(:), vrev_p(:), 'type', 'Pearson');
  ScorrV(ij) = corr(vrev(:), vrev_p(:), 'type', 'Spearman'); 
  KcorrV(ij) = corr(vrev(:), vrev_p(:), 'type', 'Kendall');
  %disp('size rxn_exp_p:');
  %disp(size(rxn_exp_p));
  %disp('size vrev_p:');
  %disp(size(vrev_p));
  PcorrEV(ij) = corr(rxn_exp_rev_p(~isnan(rxn_exp_rev_p)), abs(vrev_p(~isnan(rxn_exp_rev_p))), 'type', 'Pearson');
  ScorrEV(ij) = corr(rxn_exp_rev_p(~isnan(rxn_exp_rev_p)), abs(vrev_p(~isnan(rxn_exp_rev_p))), 'type', 'Spearman'); 
  KcorrEV(ij) = corr(rxn_exp_rev_p(~isnan(rxn_exp_rev_p)), abs(vrev_p(~isnan(rxn_exp_rev_p))), 'type', 'Kendall');
  Vactive(ij) = sum(abs(vrev_p)>Zthresh);
  vrevSign_p = signThresh(vrev_p, Zthresh);
  VdeltaSign(ij) = sum(vrevSign_p .* vrevSign < 0);
  VOnOff(ij) = sum((abs(vrevSign_p) + abs(vrevSign)) == 1);
end



