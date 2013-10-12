function testFALCON(FDBG)
% Runs some analysis and/or tests on FALCON,
% to better understand if the algorithm is 
% working correctly. 
%
% For descriptions of individual models see 
% the appropriate model-generating .m file.
%
%INPUT
% FDBG    Whether or not to turn on extra debugging info while 
%         running FALCON. This will increase the runtime.
%
% TODO: 
%
% setting equal forward/reverse rxns to zero early on can be a problem.
% explore alternatives: 
% 1) only do this if the objective value doesn't decrease
% 2) also it may only be useful at the end of the simulation anyway
%
%
% Possible alternative but non-simple algorithm to find loops? Using
% a combination of thermodynamics and inuition/heuristic -- see reviewed
% paper.
%
% Modify InTriangleOut to be a branch mode (just a couple of lines)
%
% Add in some expression tests for computeMinDisj/mindisj
%
% Test reaction groups somehow
%
% Test nvar, zvar somehow
%
% !!! Test problems with nan values. !!!
%
% Test regularization? In the triangle model, it shouldn't have an effect
% I think.

% Brandon Barker    Oct 10, 2013


EXPCON = 1;
REG = 0;

%Get directory information
thisScript = which('testFALCON');
testDir = fileparts(thisScript);
expDir = [testDir '/expression'];

disp('*** Testing FALCON with InTriangleOut model ***');

m = makeTestModel_InTriangleOut();
[mI, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(m);
[nmets, nrxns] = size(mI.S);
ngenes = length(mI.genes);
rxn_exp = nan*ones(nrxns, 1);
rxn_exp_sd = nan*ones(nrxns, 1);
rxn_rule_group = zeros(nrxns, 1);
for i = 1:ngenes
    rxn_rule_group(2*i - 1) = 2*i - 1;
    rxn_rule_group(2*i) = 2*i - 1;
    rxn_exp(2*i - 1) = 1;
    rxn_exp(2*i) = 1;
    rxn_exp_sd(2*i - 1) = 1;
    rxn_exp_sd(2*i) = 1;
end
% Last one is irreversible:
rxn_exp(2*i) = nan;
rxn_exp_sd(2*i) = nan;
for j = (2*i):nrxns
    rxn_rule_group(j) = j;
end

disp('Check that output from computeMinDisj is as expected:')
[CMD_rxn_exp, CMD_rxn_exp_sd, CMD_rxn_rule_group] = ...
    computeMinDisj(mI, [expDir '/InTriangleOut_All_1.csv'], -1, FDBG);

if all(rxn_exp(~isnan(rxn_exp)) == CMD_rxn_exp(~isnan(rxn_exp))) && ...
    all(rxn_exp(~isnan(CMD_rxn_exp)) == CMD_rxn_exp(~isnan(CMD_rxn_exp)))
    disp('Test succeeded for computeMinDisj rxn_exp');
else
    disp('Test FAILED for computeMinDisj rxn_exp');
    disp(rxn_exp)
    disp(CMD_rxn_exp)
end
if all(rxn_exp_sd(~isnan(rxn_exp_sd)) == ...
CMD_rxn_exp_sd(~isnan(rxn_exp_sd))) && ...
all(rxn_exp_sd(~isnan(CMD_rxn_exp_sd)) == ...
CMD_rxn_exp_sd(~isnan(CMD_rxn_exp_sd)))
    disp('Test succeeded for computeMinDisj rxn_exp_sd');
else
    disp('Test FAILED for computeMinDisj rxn_exp_sd');
end
if all(rxn_rule_group == CMD_rxn_rule_group)
    disp('Test succeeded for computeMinDisj rxn_rule_group');
else
    disp('Test FAILED for computeMinDisj rxn_rule_group');
end

[v_solirrev, corrval, nvar, v_all] =         ...
    falcon(mI, rxn_exp, rxn_exp_sd,          ...
           rxn_rule_group, REG, 0, EXPCON, FDBG);

v_solirrev';
v_solrev = convertIrrevFluxDistribution(m, v_solirrev, matchRev)'

if (v_solrev(1) > 0) && (v_solrev(2) > 0) && (v_solrev(7) > 0)
    disp('Test succeeded for linear pathway fluxes related to expression.');
else
    disp('Test FAILED for linear pathway fluxes related to expression.');
end

% It should be the case (unfortunately for now), that,
% given no conflicting demands, expression on a futile cycle will cause
% flux through the cycle.
if ~all(v_solrev(3:4))
    disp('Test FAILED for non-zero fluxes in a cycle with expression.');
else
    disp('Test succeeded for non-zero fluxes in a cycle with expression.');
end

disp(' ');
disp(' ');


% Same thing, but without expression on the cycle reactions
% F3 and F4.
rxn_exp_0loop = rxn_exp;
rxn_exp_0loop(5:8) = 0;

[v_solirrev_0loop, corrval, nvar, v_all] =   ...
      falcon(mI, rxn_exp_0loop, rxn_exp_sd,  ...
             rxn_rule_group, REG, 0, EXPCON, FDBG);

v_solirrev_0loop';
v_solrev_0loop = convertIrrevFluxDistribution(m, v_solirrev_0loop, matchRev)'

% It should be the case that reactions F_3 and F_4 are zero, but not F_2:
if any(v_solrev_0loop(3:4)) && ~any(v_solrev_0loop(2))
    disp('Test FAILED for zero fluxes in a 0-expression cycle.');
else
    disp('Test succeeded for zero fluxes in a 0-expression cycle.');
end

% In branch model, do a series of tests that ensure increasing
% flux with increasing expression:
% a <= b <= c <= ...

mBranch = m;
mBranch.lb(2) = 0;
mBranch.ub(3:4) = 0;
mBranch.rev(2:4) = 0;
rxn_exp_branch = rxn_exp;
rxn_exp_branch(2) = 5;
%for i = 0:10
%    rxn_exp_branch(3:4) = 
%end