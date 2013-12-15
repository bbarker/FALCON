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


% best to have an optional argument to load these all from a .mat file

% Note iIN800 and iFF708 don't appear to have valid Gene labels
yeastModFiles = {'iND750.xml', 'iMM904_flux_orig.xml', ...
                 'yeast_5.21_MCISB.xml', 'yeast_6.06_cobra.xml', ...
                 'yeast_7.00_cobra.xml'};

%include Ecoli core model
ecoliModFiles = {'ecoli_core_model.xml', 'iJR904.xml', 'iAF1260.xml' ...
                 'JO1366.xml' }

% Now save to .mat file if .mat not provided as an argument.

yeastExpFile = 
ecoliExpFile = 
humanExpFile = 



%Best to remove any directionality constraints not matching rev
%to guarantee more similarity between metrics 
%removeEnzymeIrrevs.m?