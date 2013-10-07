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
rxn_exp = zeros(1, length(nrxns));
rxn_exp_sd = ones(1, length(nrxns));
rxn_rule_group = zeros(1, length(nrxns));
for i = 1:ngenes
    rxn_rule_group(2*i - 1) = i
    rxn_rule_group(2*i) = i
end
disp('Check that output from computeMinDisj is as expected:')
[CMD_rxn_exp, CMD_rxn_exp_sd, CMD_rxn_rule_group] = ...
    computeMinDisj(mI, [expDir '/InTriangleOut_All_1.csv']);
disp(rxn_exp)
disp(CMD_rxn_exp)
disp(rxn_exp_sd)
disp(CMD_rxn_exp_sd)
disp(rxn_rule_group)
disp(CMD_rxn_rule_group)

disp(' ');
disp(' ');