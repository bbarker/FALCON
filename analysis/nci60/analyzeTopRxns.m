function [signMatch originalSign corrals originalCorral rxnsOfInt]= analyzeTopRxns (rec2, expressionFile, numTops, excFlag)
%this function looks at all possible combinations of the two
%irreversible directions for a certain number of reversible
%reactions with the highest expression levels

%add analyze_VSol...?
%INPUTS
%   rec2 is the human recon 2 model
%
%   expressionFile is the tissue being analyzed (.csv)
%
%   numTops is the number of reactions to be analyzed (up to 16 rxns)
%
%OPTIONAL INPUTS
%   excFlag is 1 if looking at top non-exchange rxns. Defaults to 0. 
%
%OUTPUTS
%       the index of any of the 1st four outputs can be used to find the combination
%       of forward (0) or backward (1) rxns by computing the vector
%       rxnsOfInt(boolean(str2numvec(dec2bin(index, length(rxnsOfInt)))))
%       i.e.:
%       0 is irreversible 0 to 1000
%       1 is irreversible -1000 to 0
%
%   signMatch is the percentage of exchange reactions that have the same
%       sign as the flux from the jain table
% 
%   originalSign is the percentage of the model without making changes
%
%   corral is the corral value returned from falcon for that combination
%
%   originalCorral is the corral of the model without making changes
%
%   rxnsOfInt is a vector containing all the reactions that were analyzed
%

% Narayanan Sadagopan  11/21/2013
% Brandon Barker       11/11/2013

if ~exist('excFlag','var')
    excFlag=0;
end

%from analyzeV_solFileOneCellLine, find tissue index in table
[cellLine metArray coretable FVAvmin FVAvmax] = readJainTable();
inputFiSplit = strsplit(expressionFile, '/');
fileName = inputFiSplit{end};
CLidx = -1;
for i = 1:length(cellLine)
    CL = convertExpressionFileName(cellLine{i});
    empty = isempty(regexp(fileName, ['^' CL]));
    if ~empty
        CLidx = i;
        break;
    end
end
if CLidx <= 0
    disp('Error: cell line not found.');
    return;
end

%declaring variables
rec2Irr = convertToIrreversible(rec2);
[rxn_exp] = computeMinDisj(rec2Irr,expressionFile);
[virrev vrev corral]=runFalcon(rec2,expressionFile,0.01,0,0);
[rxnExpSort ind] = sort(rxn_exp,'descend');
selExc = findExcRxns(rec2Irr);
rxnsOfInt = zeros(1,numTops);
count = 1;
rxnCount = 1;

%find "numTop" top reactions
if (excFlag)
    for x = 1:length(ind)
        if (rec2Irr.rev(ind(x)) && rxnExpSort(x)>0 ...
                && ~selExc(ind(rxnCount)))
            rxnsOfInt(count) = ind(x);
            count = count + 1;
            if (count>numTops)
                break;
            end
        end
    end
else
    for x = 1:length(ind)
        if (rec2Irr.rev(ind(x)) && rxnExpSort(x)>0)
            rxnsOfInt(count) = ind(x);
            count = count + 1;
            if (count>numTops)
                break;
            end
        end
    end
end

%length check and declaring size of outputs
len = length(rxnsOfInt);
if len > 16
    disp('Too many reaction sets to analyze!');
    return;
end
signMatch = zeros(1, 2^len);
corrals = zeros(1, 2^len);

%get metabolite and exchange reaction conversions
%get top 6 fluxes in jain table
excMetIds = loadJainMetsToExcIdxs(metArray, rec2);
key = keys(excMetIds);
val = values(excMetIds);
[coreSorted coreInd]= sort(abs(coretable(:,CLidx)),'descend');
keyCount = 1;
for y = 1 : length(coreInd)
    for z = 1 : length(key)
        if (strcmp(metArray(coreInd(y)),key(z)))
            key2{keyCount} = key(z);
            val2{keyCount} = val{z};
            keyCount = keyCount + 1;
        end
    end
    if (keyCount>6)
        break;
    end
end
            
%find originalSign
signCount = 0; 
for y = 1 : length(metArray)
    for z = 1 : length(key2)
        if (strcmp(key2{z}, metArray(y)))
            if (length(val2{z})==1)
                if ((coretable(y,CLidx) > 0 && vrev(val2{z}) > 0) ...
                        ||(coretable(y,CLidx) < 0 && vrev(val2{z}) < 0)...
                        ||(coretable(y,CLidx) == 0 && vrev(val2{z}) == 0))
                    signCount = signCount + 1;
                end
                break;
            else
                temp = 0;
                for i = 1:length(val2{z})
                    temp = temp + vrev(val2{z}(i));
                end
                if ((coretable(y,CLidx) > 0 && temp > 0) ...
                        ||(coretable(y,CLidx) < 0 && temp < 0)...
                        ||(coretable(y,CLidx) == 0 && temp == 0))
                    signCount = signCount + 1;
                end
                break;
            end
        end
    end
end
originalSign = signCount/length(key2);
originalCorral = corral;

%run for all combinations of 0 to 1000 and -1000 to 0
parfor x = 1 : 2^len
    recMod = rec2;
    aX = boolean(str2numvec(dec2bin(x - 1, len)));
    recMod.rev(rxnsOfInt) = 0;
    recMod.lb(rxnsOfInt(~aX)) = 0;
    recMod.ub(rxnsOfInt(~aX)) = 1000;
    recMod.lbrxnsOfInt(rxnsOfInt(aX)) = -1000;
    recMod.ub(rxnsOfInt(aX)) = 0;
    [virrev vrev corral] = runFalcon(recMod, expressionFile, 0.01, false, 0);
    signCount = 0; 
    for y = 1 : length(metArray)
        for z = 1 : length(key2)
            if (strcmp(key2{z}, metArray(y)))
                if (length(val2{z})==1)
                    if ((coretable(y,CLidx) > 0 && vrev(val2{z}) > 0) ...
                            ||(coretable(y,CLidx) < 0 && vrev(val2{z}) < 0)...
                            ||(coretable(y,CLidx) == 0 && vrev(val2{z}) == 0))
                        signCount = signCount + 1;
                    end
                    break;
                else
                    temp = 0;
                    for i = 1:length(val2{z})
                        temp = temp + vrev(val2{z}(i));
                    end
                    if ((coretable(y,CLidx) > 0 && temp > 0) ...
                            ||(coretable(y,CLidx) < 0 && temp < 0)...
                            ||(coretable(y,CLidx) == 0 && temp == 0))
                        signCount = signCount + 1;
                    end
                    break;
                end
            end
        end
    end
    signMatch(x) = signCount/length(key2);
    corrals(x) = corral;
end