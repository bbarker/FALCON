function meanTime = 
    methodTimer(model, genedata_filename, method, gene_to_scale, flux_to_scale, nReps)
%
% meanTime is calculated as average falcon or lee run + one instance of complexation
%

expCon = false;
minFit = 0.0;
regC = 0;

nrxns = length(model.rxns);
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

%
%
if strcmp(method, 'FALCON')
%
%

mdStart = tic;
[rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
    computeMinDisj(modelIrrev, genedata_filename);
mdTime = toc(mdStart);

falconStart = tic;
falconMulti(modelIrrev, nReps, rxn_exp_md,              ...
            rxn_exp_sd_md, rxn_rule_group, 'rc', regC,  ...
            'minFit', minFit, 'EXPCON', expCon);
meanTime = toc(falconStart)/nReps + mdTime;

end % of if strcmp(method, 'FALCON')

%
%
if strcmp(method, 'eMoMA')
%
%

expStart = tic;
[r_lee, rs_lee, r_miss] = geneToRxn(model, genedata_filename);
expTime = toc(expStart);

leeStart = tic;
dataToFluxFixMulti(m, r_lee, rs_lee, nReps);
meanTime = toc(leeStart)/nReps + expTime;


end % if strcmp(method, 'eMoMA')
