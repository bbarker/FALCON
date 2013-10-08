function testFALCON()
% Runs some analysis and/or tests on FALCON,
% to better understand if the algorithm is 
% working correctly.
%
% For descriptions of individual models see 
% the appropriate model-generating .m file.


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
for j = (2*i + 1):nrxns
    rxn_rule_group(j) = j;
end
disp('Check that output from computeMinDisj is as expected:')
[CMD_rxn_exp, CMD_rxn_exp_sd, CMD_rxn_rule_group] = ...
    computeMinDisj(mI, [expDir '/InTriangleOut_All_1.csv']);

if all(rxn_exp(~isnan(rxn_exp)) == CMD_rxn_exp(~isnan(rxn_exp))) && ...
   all(rxn_exp(~isnan(CMD_rxn_exp)) == CMD_rxn_exp(~isnan(CMD_rxn_exp)))
    disp('Test succeeded for computeMinDisj rxn_exp');
else
    disp('Test FAILED for computeMinDisj rxn_exp');
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

disp(' ');
disp(' ');