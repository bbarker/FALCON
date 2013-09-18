function analyzeMultiPerturbation(model, expFileDir)
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
%OUTPUT      a file in 'outputDir' (defined below) for each
%            cell line.
%
% Brandon Barker 09/16/13
%

%Get a list of subdirectories in the specified directory
pertPaths = setdiff(strsplit(genpath(expFileDir),':'), {expFileDir});
pertPaths = pertPaths(boolean(cellfun(@length, pertPaths)));

%expressionFile = convertExpressionFileName(CL);

statsMegaCell = cell([1, length(pertPaths)]);
maxStasLen = 0;
numCellLines = 0;
for i = 1:length(pertPaths)
    fluxFileList = dir([pertPaths{i} '/*out']);
    fluxFileList = struct2cell(fluxFileList);
    fluxFileList = fluxFileList(1,:);
    expSubDir = pertPaths{i};
    simParams = directoryLabelParse(expSubDir,'_');
    numCellLines = length(fluxFileList);
    statsCell = cell([1, numCellLines]);
    parfor j = 1:numCellLines
        fluxFile = [pertPaths{i} '/' fluxFileList{j}];
	statsCell{j} = analyzeV_solFileOneCellLine(fluxFile);
    end
    statsMegaCell{i} = statsCell;
end


flatAnalysisFI = fopen([expFileDir 'flatAnalysis.csv'], 'w');

for i = 1:numCellLines
    for j = 1:length(pertPaths)
        statsCell = statsMegaCell{j};
	statsLen = length(keys(statsCell{i}));
	disp(['stats length: ' num2str(statsLen)]);
        valString = '';

        for k = 1:(statsLen - 3) % 3 fields are non-integer
	    disp(num2str(k));
	    disp(keys(statsCell{i}));
	    statsArray = statsCell{i}(num2str(k));
            valString = [valString '\t' num2str(statsArray(1))];
        end
        %fprintf(flatAnalysisFI,,A1,...,An)
    end
end

fclose(flatAnalysisFI);
