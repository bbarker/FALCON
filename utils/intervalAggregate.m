function [x, y] = intervalAggregate(Xdata, Ydata, aggFun, intStep, intOverlap)
% intOverlap in [0, 1); 0 for no overlap of intervals, etc.
% intStep    this is the size of the interval being aggregated.

minX = min(Xdata);
maxX = max(Xdata);

minY = min(Ydata);
maxY = max(Ydata);

intInc = intOverlap*intStep; %How far we advance each iteraction.
nInt = ceil((maxX-minX)/intInc); %Number of aggregations

parfor i = 1:nInt
    xStart = minX + (i-1)*intInc;
    xEnd   = xStart + intStep; 
    intervalIndices = find((Xdata >= xStart) & (Xdata <= xEnd));
    x(i) = aggFun(Xdata(intervalIndices));
    y(i) = aggFun(Ydata(intervalIndices));
end

