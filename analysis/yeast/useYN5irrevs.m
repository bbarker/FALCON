function ynOut = useYN5irrevs(yn5, ynNew)
% since we require some irreversible reactions to start ...
% get these from the published paper.

ynOut = ynNew;
for i = 1:length(yn5.rxns)
    if (yn5.rev(i) == 0) && numel(strfind(yn5.rxnNames{i}, 'isa')) == 0
        rxn = yn5.rxnNames{i};
        idxIn7 = find(strcmp(ynNew.rxnNames, rxn));
        if idxIn7
            ynOut.rev(idxIn7) = 0;
            disp(ynOut.rxnNames(idxIn7));
            if yn5.lb(i) == 0
                ynOut.lb(idxIn7) = 0;
            end
            if yn5.ub(i) == 0
                ynOut.ub(idxIn7) = 0;
            end
        else
            disp(['Could not find ' rxn]);
        end
    end
end

