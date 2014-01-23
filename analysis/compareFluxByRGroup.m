function [v_md, v_md_nogrp, nnanDiffFlux] = ...
    compareFluxByRGroup(model, expFile, nReps, fLabel)

%Simple script to compare fluxes output from the FALCON algorithm
%using data on reactions that share enzyme complex or not.
%
%OUTPUT
%
%nnanDiffFlux        Number of differences in FALCON flux for
%
%v_grp               Flux esitmated with group data (default)
% 
%v_no_grp            Flux estimate without reaction group data.           
%
%These are all saved to a .mat file (see save() below).

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

nRxns = length(model.rxns);
nRxnsIrr = length(modelIrrev.rxns);

[r_mdIrr, rs_mdIrr, r_groupIrr] = computeMinDisj(modelIrrev, expFile);
fake_r_groupIrr = [1:nRxnsIrr]';

[v_mdIrr, ~, ~, ~, ~, ~, v_mdIrr_s] =  ...
    falconMulti(modelIrrev, nReps, r_mdIrr, rs_mdIrr, r_groupIrr);

[v_md_nogrpIrr, ~, ~, ~, ~, ~, v_md_nogrpIrr_s] =  ...
    falconMulti(modelIrrev, nReps, r_mdIrr, rs_mdIrr, fake_r_groupIrr);


v_md = convertIrrevFluxDistribution(v_mdIrr, matchRev);
v_md_nogrp = convertIrrevFluxDistribution(v_md_nogrpIrr, matchRev);
v_md_s = abs(convertIrrevFluxDistribution(v_mdIrr_s, matchRev));
v_md_nogrp_s = abs(convertIrrevFluxDistribution(v_md_nogrpIrr_s, matchRev));

nnanDiffFlux = sum(isDiffExp(v_md, v_md_nogrp, false(nRxns, 1))); 

save(['FluxByGroupComp_' model.description expFile fLabel '_' num2str(nReps) ...
      '.mat'], 'v_md', 'v_md_nogrp', 'v_md_s', 'v_md_nogrp_s', 'nnanDiffFlux');

function d = isDiffExp(r1, r2, rmg)
d = columnVector(boolean(abs(r1 - r2) > 1e-4) & (~rmg))';
% end of isDiffExp
