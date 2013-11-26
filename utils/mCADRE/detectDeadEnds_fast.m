function [deadEnd,onlyProd,onlyCons,bothProdCons] = detectDeadEnds_fast(model)
%detectDeadEnds Detect dead ends in the model
%
% [deadEnd,onlyProd,onlyCons,bothProdCons] = detectDeadEnds(model)
%
% Markus Herrgard

deadEnd = [];
onlyProd = [];
onlyCons = [];
bothProdCons = [];
for metNo = 1:length(model.mets)
    bothRxns = find(model.S(metNo,:) ~= 0 & model.rev(:)');
    prodRxns = union(find(model.S(metNo,:) > 0),bothRxns);
    consRxns = union(find(model.S(metNo,:) < 0),bothRxns);
    if (isempty(consRxns))
        deadEnd(end+1) = metNo;
        onlyProd(end+1) = metNo;
    elseif (isempty(prodRxns))
        deadEnd(end+1) = metNo;
        onlyCons(end+1) = metNo;
    else
        if (length(prodRxns) == 1 & length(consRxns)==1 & prodRxns == consRxns)
            deadEnd(end+1) = metNo;
            bothProdCons(end+1) = metNo;
        end
    end 
end
