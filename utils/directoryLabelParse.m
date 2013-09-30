function params = directoryLabelParse(aPath, delim, endDelim)
%
% Returns a vector of parameters (or perturbations) where the
% parameters are numeric values separated by delim in the
% leaf subdirectory.
%
%INPUT
% aPath   a relative or absolute path to a folder that is all
%         numeric values separated by delim
%
% delim   a non-numeric string delimiter
%
%OUTPUT
% params vector of numeric values
%

folders = strsplit(aPath, '/');
subDir = folders{end};
if exist('endDelim', 'var')
    subDirP = strsplit(subDir, endDelim);
    subDir = subDirP{1};
end
params = strsplit(subDir, delim);
params = cellfun(@str2num, params);