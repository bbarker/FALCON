function [scores rxnsOfInt signMatch]= analyzeTopRxns (rec2, expressionFile, numTops, excFlag)
%this function looks at all possible combinations of the two
%irreversible directions for a certain number of reversible
%reactions with the highest fluxes

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
%   scores is the absolute value of the differences in fluxes between
%       the jain table and the falcon generated fluxes. 
%       the index in scores can be used to find the combination
%       of forward (0) or backward (1) rxns by computing the vector
%       rxnsOfInt(boolean(str2numvec(dec2bin(index, length(rxnsOfInt)))))
%       i.e.:
%       0 is irreversible 0 to 1000
%       1 is irreversible -1000 to 0
%   rxnsOfInt is a vector containing all the reactions that were analyzed
%

% Narayanan Sadagopan  11/10/2013
% Brandon Barker       11/11/2013

if ~exist('excFlag','var')
    excFlag=0;
end

%from analyzeV_solFileOneCellLine
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

[virrev vrev] = runFalcon(rec2, expressionFile, 0.01, false, 0);
[vrevSort ind] = sort(abs(vrev),'descend');
selExc=findExcRxns(rec2);

rxnsOfInt = zeros(1,numTops);
count = 1;
rxnCount = 1;
%selects "numTops" top reactions
if (excFlag)
    while (count <= numTops)
        if (rec2.rev(ind(rxnCount)) && ~selExc(ind(rxnCount)))
            rxnsOfInt(count) = ind(rxnCount);
            count = count + 1;
        end
        rxnCount = rxnCount + 1;
    end
else
    while (count <= numTops)
        if (rec2.rev(ind(rxnCount)))
            rxnsOfInt(count) = ind(rxnCount);
            count = count + 1;
        end
        rxnCount = rxnCount + 1;
    end
end

len = length(rxnsOfInt);
if len > 16
    disp('Too many reaction sets to analyze!');
    return;
end
scores = zeros(1, 2^len);
signMatch = zeros(1, 2^len);

excMetIds = loadJainMetsToExcIdxs(metArray, rec2);
key = keys(excMetIds);
val = values(excMetIds);

parfor x = 1 : 2^len
    recMod = rec2;
    aX = boolean(str2numvec(dec2bin(x - 1, len)));
    recMod.rev(rxnsOfInt) = 0;
    recMod.lb(rxnsOfInt(~aX)) = 0;
    recMod.ub(rxnsOfInt(~aX)) = 1000;
    recMod.lbrxnsOfInt(rxnsOfInt(aX)) = -1000;
    recMod.ub(rxnsOfInt(aX)) = 0;
    [virrev vrev] = runFalcon(recMod, expressionFile, 0.01, false, 0);
    tempScore = 0;
    signCount = 0; 
    for y = 1 : length(metArray)
        for z = 1 : length(key)
            if (strcmp(key(z), metArray(y)))
                if (length(val{z})==1)
                    tempScore = tempScore + abs(coretable(y,CLidx)-vrev(val{z}));
                    if ((coretable(y,CLidx) > 0 && vrev(val{z}) > 0) ...
                            ||(coretable(y,CLidx) < 0 && vrev(val{z}) < 0))
                        signCount = signCount + 1;
                    end
                    break;
                else
                    temp = 0;
                    for i = 1:length(val{z})
                        temp = temp + vrev(val{z}(i));
                    end
                    tempScore = tempScore + abs(coretable(y,CLidx)-temp);
                    if ((coretable(y,CLidx) > 0 && temp > 0) ...
                            ||(coretable(y,CLidx) < 0 && temp < 0))
                        signCount = signCount + 1;
                    end
                    break;
                end
            end
        end
    end
    signMatch(x) = signCount/length(key);
    scores(x) = tempScore;
end
