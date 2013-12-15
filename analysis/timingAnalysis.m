function [timeMat iterMat] = timingAnalysis(nReps, models)
%
% Set this to the location of the SBML models directory. 
myModelDir = '/home/brandon/FBA/models';
%
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
if ~exist('models', 'var')
    origDir = pwd();
    cd(myModelDir);
    models = {};
    % Note iIN800 and iFF708 don't appear to have valid Gene labels
    yeastModFiles = {'iND750', 'iMM904_flux_orig', 'yeast_5.21_MCISB', ...
                     'yeast_6.06_cobra', 'yeast_7.00_cobra'};

    for i = 1:length(yeastModFiles)
        models{end+1} = readCbModel([yeastModFiles{i} '.xml']);
        models{end} = removeEnzymeIrrevs(models{end});
    end
    %include Ecoli core model
    ecoliModFiles = {'ecoli_core_model', 'iJR904', 'iAF1260', 'iJO1366'};
    for i = 1:length(ecoliModFiles)
        models{end+1} = readCbModel([ecoliModFiles{i} '.xml']);
        models{end} = removeEnzymeIrrevs(models{end});
    end
    % ? include recon 2.1 if on time?
    humanModFiles = {'recon1', 'recon2model.v02'};
    for i = 1:length(humanModFiles)
        models{end+1} = readCbModel([humanModFiles{i} '.xml']);
        models{end} = initializeRecon2(models{end});
        %removeEnzymeIrrevs handled in initializeRecon2
    end
    %We can use a set of tissue specific models for human recon 1 and 2
    %but this might be overkill for the original implementation,
    %and may introduce other complications from tissue-specific model
    %creation

    % Now save to .mat file if .mat not provided as an argument.
    cd(origDir);
    save('timingModels.mat', 'models');
end


%yeastExpFile = 
%ecoliExpFile = 
%humanExpFile = 



%Best to remove any directionality constraints not matching rev
%to guarantee more similarity between metrics 
%removeEnzymeIrrevs.m?