function meanTime = ...
    methodTimer(model, genedata_filename, method, nReps, gene_to_scale)
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

[~, ~, ~, ~, fTime] = falconMulti(modelIrrev, nReps, rxn_exp_md, ...
            rxn_exp_sd_md, rxn_rule_group, 'rc', regC,             ...
            'minFit', minFit, 'EXPCON', expCon);

meanTime = fTime + mdTime;

end % of if strcmp(method, 'FALCON')

%
%
if strcmp(method, 'eMoMA')
%
%

m = model;

expStart = tic;
[r_lee, rs_lee, r_miss] = geneToRxn(m, genedata_filename);
%
uptake          = find(strcmp(gene_to_scale, m.rxnNames));
rs_lee          = rs_lee/r_lee(uptake);
r_lee           = r_lee/r_lee(uptake);
%
m.lb(uptake)	= 1;
m.ub(uptake)	= 1;
%
expTime = toc(expStart);

leeStart = tic;
[~, ~, ~, ~, lTime] = ...
dataToFluxFixMulti(m, r_lee, rs_lee, nReps);
meanTime = lTime + expTime;

end % if strcmp(method, 'eMoMA')
