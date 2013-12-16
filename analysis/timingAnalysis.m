function [timeMat iterMat] = timingAnalysis(nReps, method, models)
%
% Set this to the location of the SBML models directory. 
myModelDir = '/home/brandon/FBA/models';
% Set these to the expression files:
yeastExpFile = '/home/brandon/FBA/models/Analysis/AdaptiveMOMA/LeeExpFlux/genedata_75.txt';
ecoliExpFile = '/home/brandon/FBA/models/Analysis/AdaptiveMOMA/LeeExpFlux/ecoli_MOPSCompGluc.expSTD.csv';
humanExpFile = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/nci60rseq_thresh/0_0.0/K562.csv';

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


yeastModFiles = {'iND750', 'iMM904_flux_orig', 'yeast_5.21_MCISB', ...
                 'yeast_6.06_cobra', 'yeast_7.00_cobra'};
ecoliModFiles = {'ecoli_core_model', 'iJR904', 'iAF1260', 'iJO1366'};
humanModFiles = {'recon1', 'recon2model.v02'};
%Note that 2.02 has some bugs, we replace it with 2.03, only available
%as a .mat file on humanmetabolism.org.
modFileNames = cellfun(@(x) [x '.xml'], {yeastModFiles{:}, ecoliModFiles{:}, ...
                  humanModFiles{:}}, 'UniformOutput', false);
nModels = length(modFileNames);

nYmod = length(yeastModFiles);
nECmod = length(ecoliModFiles);
nHmod = length(humanModFiles);

timeMat = nan*ones(nModels, nReps);
iterMat = nan*ones(nModels, nReps); 

expFileNames(1 : nYmod) = deal({yeastExpFile});
expFileNames(nYmod + 1 : nYmod + nECmod) = deal({ecoliExpFile});
expFileNames(nYmod + nECmod + 1 : nYmod + nECmod + nHmod) = deal({humanExpFile});


% best to have an optional argument to load these all from a .mat file


if ~exist('models', 'var')
    origDir = pwd();
    cd(myModelDir);
    models = {};
    % Note iIN800 and iFF708 don't appear to have valid Gene labels

    for i = 1:nYmod
        models{end+1} = readCbModel([yeastModFiles{i} '.xml']);
        models{end} = removeEnzymeIrrevs(models{end});
    end
    %include Ecoli core model

    for i = 1:nECmod
        models{end+1} = readCbModel([ecoliModFiles{i} '.xml']);
        models{end} = removeEnzymeIrrevs(models{end});
    end
    % ? include recon 2.1 if on time?
    for i = 1:nHmod
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

for i = 1:nModels
    mI = convertToIrreversible(models{i});
    iVec = nan*ones(1, nReps);
    tVec = nan*ones(1, nReps);
    if strcmp(method, 'FALCON')
        parfor j = 1:nReps
            t0 = tic();
            [r, r_s, r_g] = computeMinDisj(mI, expFileNames{i});
            [~, ~, ~, ~, ~, fIter] = falcon(mI, r, r_s, r_g);
            tVec(j) = toc(t0);
            iVec(j) = fIter;
        end
    elseif strcmp(method, 'Lee')
        parfor j = 1:nReps
            t0 = tic();
            [r, r_s, r_misG] = geneToRxn(models{i}, expFileNames{i});
            [~, lIter] = dataToFluxFix(models{i}, r, r_s);
            tVec(j) = toc(t0);
            iVec(j) = lIter;
        end 
    end
    timeMat(i, :) = tVec;
    iterMat(i, :) = iVec;
end

timeMat = timeMat';
iterMat = iterMat';

%Best to remove any directionality constraints not matching rev
%to guarantee more similarity between metrics 
%removeEnzymeIrrevs.m?