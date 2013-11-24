function boundsRev = getBoundsRev(m)
% A fast function to determine reversibility in an irreversible
% model based on mathematical model information:
% m.S and m.ub

nrxns = length(m.rxns);
nmets = length(m.mets);
boundsRev = zeros(nrxns, 1);

for i = 1 : (nrxns - 1)
    if all(m.S(:, i) ==  -m.S(:, i+1))
        if (m.ub(i) > 0) && (m.ub(i+1) > 0)
            boundsRev(i : (i+1)) = 1;
            i = i + 1;
        end
    end
end
