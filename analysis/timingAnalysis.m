function [timeMat iterMat] = timingAnalysis( nReps)

% Takes a cell array of models, and uses
% the internal lookup table to find specified
% data files.

% In the future, if the model doesn't have expression data available, use
% some form of ortholog or EC reaction search to sample from available
% data. 

% For now, we might be able to use reaction lookup name (assuming 
% standard reaction ids or names, and if missing, using a random
% expression value or median expression value.

% For the plot, we can most likely plot the points we have
% real data for in a different color.