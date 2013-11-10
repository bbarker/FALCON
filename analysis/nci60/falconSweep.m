function [scores a]= falconSweep (rec2,tissueName,rxnsOfInt)
%this function looks at all possible combinations of 
%reversible and irreversible for the input reactions

%INPUTS
%   rec2 is the human recon 2 model
%   tissueName is the tissue being analyzed (.csv)
%   rxnsOfInt is a vector containing the reactions to be analyzed (up
%       to 6 rxns)
%
%OUTPUTS
%   a is the list of all combinations run
%   scores is whether glucose is transported in or out

len = length(rxnsOfInt);

v = [0 1 2];
%0 is irreversible 0 to 1000
%1 is irreversible -1000 to 0
%2 is reversible -1000 to 1000

if (len==1)
    a=allcomb(v);
elseif (len==2)
    a=allcomb(v,v);
elseif (len==3)
    a=allcomb(v,v,v);
elseif (len==4)
    a=allcomb(v,v,v,v);
elseif (len==5)
    a=allcomb(v,v,v,v,v);
else
    a=allcomb(v,v,v,v,v,v);
end

scores = zeros(1,length(a(:,1)));

for x = 1:length(a(:,1))
    recMod = rec2;
    for y = 1:length(a(x,:))
        if (a(x,y)==0)
            recMod.lb(rxnsOfInt(y)) = 0;
            recMod.ub(rxnsOfInt(y)) = 1000;
            recMod.rev(rxnsOfInt(y)) = 0;
        elseif (a(x,y)==1)
            recMod.lb(rxnsOfInt(y)) = -1000;
            recMod.ub(rxnsOfInt(y)) = 0;
            recMod.rev(rxnsOfInt(y)) = 0;
        else
            recMod.lb(rxnsOfInt(y)) = -1000;
            recMod.ub(rxnsOfInt(y)) = 1000;
            recMod.rev(rxnsOfInt(y)) = 1;
        end
    end
    [virrev vrev] = runFalcon(recMod,tissueName,0.01,true,0);
    if (vrev(1388)<0)
        scores(x) = 1;
    elseif (vrev(1388)>0)
        scores(x) = 2;
    end
end
