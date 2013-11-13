function [scores]= falconSweep (rec2, expressionFile, rxnsOfInt, targetRxn, compZero)
%this function looks at all possible combinations of the two
%irreversible directions for the input reactions

%INPUTS
%   rec2 is the human recon 2 model
%
%   expressionFile is the tissue being analyzed (.csv)
%
%   rxnsOfInt is a vector containing the reactions to be analyzed (up
%       to 16 rxns)
%
%   targetRxn is the reaction to who's flux is being analyzed
%
%   compZero is 0 if looking for less than zero fluxes, 1 if looking for 
%       greater than zero fluxes
%
%OUTPUTS
%   scores is whether glucose is transported in or out
%       the index in scores can be used to find the combination
%       of forward (0) or backward (1) rxns by computing the vector
%       rxnsOfInt(boolean(str2numvec(dec2bin(index, length(rxnsOfInt)))))
%       i.e.:
%       0 is irreversible 0 to 1000
%       1 is irreversible -1000 to 0

% Narayanan Sadagopan  11/10/2013
% Brandon Barker       11/11/2013

len = length(rxnsOfInt);

if len > 16
    disp('Too many reaction sets to analyze!');
    return;
end

scores = zeros(1, 2^len);

parfor x = 1 : 2^len
    recMod = rec2;
    aX = boolean(str2numvec(dec2bin(x - 1, len)));
    recMod.rev(rxnsOfInt) = 0;
    recMod.lb(rxnsOfInt(~aX)) = 0;
    recMod.ub(rxnsOfInt(~aX)) = 1000;
    recMod.lbrxnsOfInt(rxnsOfInt(aX)) = -1000;
    recMod.ub(rxnsOfInt(aX)) = 0;
    [virrev vrev] = runFalcon(recMod, expressionFile, 0.01, false, 0);
    if (compZero==0)
        if (vrev(targetRxn)<0)
            scores(x)=1;
        end
    else
        if (vrev(targetRxn)>0)
            scores(x) = 1;
        end
    end
end