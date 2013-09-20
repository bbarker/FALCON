function runMultiPerturbtion(model, expFileDir, CL, envConstrain, addLabel)
%
% Calls runComparisonScript for FALCON for a single cell line,
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
% CL            If nonempty or exists, should be the name
%               of a single cell-line to run.
%
%
%OPTIONAL INPUTS
% envConstrain   'medium', 'core', or 'core_med':
%                whether or not to constraints that are 
%                medium-based CoRe-based, or both.
%
% addLabel       Label for any changes made to the model 
%                to be used in naming the output directory.
%
%OUTPUT      a file in 'outputDir' (defined below) for each
%            cell line.
%
% Brandon Barker 09/15/13
%

%Get a list of subdirectories in the specified directory
pertPaths = setdiff(strsplit(genpath(expFileDir),':'), {expFileDir});
pertPaths = pertPaths(boolean(cellfun(@length, pertPaths)));
expressionFile = convertExpressionFileName(CL);

modLabel = '';
if exist('addLabel', 'var')
    if length(addLabel) > 0
	modLabel = addLabel;
    end
end

parfor i = 1:length(pertPaths)
    expSubDir = pertPaths{i};
    runComparisonScript(model, 'FALCON', expSubDir, envConstrain, CL, modLabel);
end