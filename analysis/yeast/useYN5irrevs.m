function [ynOut, ynBoundChanges] = useYN5irrevs(yn5, ynNew, leaveOut, noDir)
% since we require some irreversible reactions to start ...
% get these from the published paper.

% Just in case it hasn't been run:
if ~exist('noDir', 'var')
    yn5   = removeEnzymeIrrevs(yn5);
    ynNew = removeEnzymeIrrevs(ynNew);
end

ynOut = ynNew;
ynLBOldIdx = [];
ynUBOldIdx = [];
ynLBOldVal = [];
ynUBOldVal = [];
idxIn7Prev = [];
if exist('leaveOut', 'var')
    idxIn7Prev = leaveOut;
end
for i = 1:length(yn5.rxns)
    if (yn5.rev(i) == 0) && numel(strfind(yn5.rxnNames{i}, 'isa')) == 0
        rxn = yn5.rxnNames{i};
        idxIn7 = find(strcmp(ynNew.rxnNames, rxn));
        idxIn7 = setdiff(idxIn7, idxIn7Prev);
        if idxIn7
            idxIn7Prev = union(idxIn7Prev, idxIn7);
            ynOut.rev(idxIn7) = 0;
            %disp(ynOut.rxnNames(idxIn7));
            if yn5.lb(i) == 0
                for i7i = 1:length(idxIn7)
                    i7 = idxIn7(i7i);
                    if ynNew.lb(i7) ~= 0
                        ynOut.lb(i7) = 0;
                        ynLBOldIdx = [ynLBOldIdx i7];
                        ynLBOldVal = [ynLBOldVal ynNew.lb(i7)];
                    end
                end
            end
            if yn5.ub(i) == 0
                for i7i = 1:length(idxIn7)
                    i7 = idxIn7(i7i);
                    if ynNew.ub(i7) ~= 0
                        ynOut.ub(i7) = 0;
                        ynUBOldIdx = [ynUBOldIdx i7];
                        ynUBOldVal = [ynUBOldVal ynNew.ub(i7)];
                    end
                end
            end
        else
            1;
            %disp(['Could not find ' rxn]);
        end
    end
end

ynBoundChanges.LBOldIdx = ynLBOldIdx;
ynBoundChanges.UBOldIdx = ynUBOldIdx;
ynBoundChanges.LBOldVal = ynLBOldVal;
ynBoundChanges.UBOldVal = ynUBOldVal;

