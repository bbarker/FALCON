function analyzeMultiPerturbtion(model, expFileDir)
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

parfor i = 1:length(pertPaths)
    expSubDir = pertPaths{i};
    simParams = directoryLabelParse(expSubDir,'_');

end

