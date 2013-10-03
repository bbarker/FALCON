function analyzeCLPermutations(model, expFileDir, colOrder, paramRange)
%
% This script should be called after careful inspection of output
% from e.g. analyzeMultiPerturbation.m and 
% analyzeMultiPerturbation_ErrorBars.m, because we need to run
% such scripts to select a good parameter ***range***. 
%
% Note: model not used yet
%
% Analyzes FALCON output from runMultiPerturbtion for a single cell line,
% and searches subdirectories for files beloging to this cell line,
% which should have directory labels corresponding to their perturbation
% and condition. If multiple perturbations apply, then a "_" should 
% be the delimiter for these perturbations in the directory name.
%
%INPUT
% model   (reversible form; the following fields are required)
%   S            Stoichiometric matrix
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         reaction ids
%   
%   The following fields for model may be required depending
%   on the method:
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%
% 
% expFileDir    Relative path to top level directory containing
%               cell-line expression files. These files are
%               tab-delimited file with a head for the columns:
%               gene (entrez gene id), mean (expression value,
%               and standard deviation (of expression).
%
%OPTIONAL INPUTS
% colOrder   vector used as sort order for 'sortrows',
%            used in determining the traverasl over the set of
%            parameters. 
% 
%OUTPUT      a file in 'outputDir' (defined below) for each
%            cell line.
%
% Brandon Barker 09/16/13
%

tab = sprintf('\t');

%Get a list of subdirectories in the specified directory
pertPaths = setdiff(strsplit(genpath(expFileDir),':'), {expFileDir});
pertPaths = pertPaths(boolean(cellfun(@length, pertPaths)));

%expressionFile = convertExpressionFileName(CL);

statsMegaCell = cell([1, length(pertPaths)]);
maxStasLen = 0;
numCellLines = 0;
fluxFileList = {};

nParams = length(directoryLabelParse(pertPaths{1}, '_', '~'));
paramMat = zeros(length(pertPaths), nParams);
for i = 1:length(pertPaths)
    expSubDir = pertPaths{i};
    paramMat(i,:) = directoryLabelParse(expSubDir, '_', '~');
end

if exist('colOrder', 'var')
    [paramsSorted, paramsIdx] = sortrows(paramMat, colOrder);
else
    [paramsSorted, paramsIdx] = sortrows(paramMat);
end

for i = 1:length(pertPaths)
    ppIdx = paramsIdx(i);
    fluxFileList = dir([pertPaths{ppIdx} '/*out']);
    fluxFileList = struct2cell(fluxFileList);
    % This assumes all the subdirectories have
    % the same file names.
    fluxFileList = sort(fluxFileList(1,:));
    numCellLines = length(fluxFileList);
    subsetsToStatsCell = cell([1, numCellLines]);
    parfor j = 1:numCellLines
        fluxFile = [pertPaths{ppIdx} '/' fluxFileList{j}];
	subsetsToStatsCell{j} = analyzeV_solFileOneCellLine(fluxFile);
    end
    statsMegaCell{ppIdx} = subsetsToStatsCell;
end


% How many ranked exo metabolites to consider
metDepth = [2 3 4 5 6 10 15]; 

% Which stats to get and their labels.
statLabel = {'Pearson', 'Spearman', 'Kendall', 'CosineSim', ...
             'L1Dist', 'UptakeSensitivity', 'ReleaseSensitivity'};

% Having this index-based label mechanism is quite unfortunate; need to think
% of how to make it intrinsic to analyzeV_solFileOneCellLine.m
% Preferably use a second container.Map for the array of values.
% allowing for intrinsic names?

dirOut = [rstrtrim(expFileDir, '/') '_stats'];
mkdir(dirOut);

for i = 1:numCellLines
    for k = 1:length(metDepth)
	dstr = num2str(metDepth(k));
        for s = 1:length(statLabel)
	    skey = statLabel{s};
            valString = '';
	    for j = 1:length(pertPaths)
                ppIdx = paramsIdx(j);
		subsetsToStatsCell = statsMegaCell{ppIdx};
		statsCell = subsetsToStatsCell{i}(dstr);
		% valString = [valString '\t' num2str(statsCell(skey))];
                valString = sprintf('%s\t%s', valString, ...
                                    num2str(statsCell(skey)));
	    end
	    valString = lstrtrim(valString, tab);
            % Probably need to use sprintf here.
	    fileName = [dirOut '/' fluxFileList{i} '_' ... 
                        statLabel{s} dstr '.csv'];
            flatAnalysisFI = fopen(fileName, 'w');
	    fprintf(flatAnalysisFI, '%s\n', valString);
            fclose(flatAnalysisFI);
        end
    end
end

