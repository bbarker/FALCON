function [x, y] = intervalAggregate(Xdata, Ydata, aggFun, intStep, intOverlap)
% intOverlap in [0, 1); 0 for no overlap of intervals, etc.
% intStep    this is the size of the interval being aggregated.

% Brandon Barker 01/21/2014
% For possible updates, see http://stackoverflow.com/questions/12556491/how-to-fit-a-curve-by-a-series-of-segmented-lines-in-matlab/21271603#21271603

minX = min(Xdata);
maxX = max(Xdata);

minY = min(Ydata);
maxY = max(Ydata);

intInc = intOverlap*intStep; %How far we advance each iteraction.
if intOverlap <= 0
  intInc = intStep;
end
nInt = ceil((maxX-minX)/intInc); %Number of aggregations

parfor i = 1:nInt
    xStart = minX + (i-1)*intInc;
    xEnd   = xStart + intStep; 
    intervalIndices = find((Xdata >= xStart) & (Xdata <= xEnd));
    x(i) = aggFun(Xdata(intervalIndices));
    y(i) = aggFun(Ydata(intervalIndices));
end

