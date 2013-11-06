function runAnalyzeFirst8Deep(model, paramDir, analysisOnly, allCL, EXPCON)
% A simple script to run 8 of 9 (since 8 is easy to do in parallel)
% of the cell lines for which we have "deep protoeomic data". Can
% optionally run and analyze all cell lines.
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
%OPTIONAL INPUTS
%
% analysisOnly     Assumes simulations of interest are finished.
%                  and will redo analysis (useful if analysis 
%                  pipeline changes at all).
%
% allCL            If true, will run all available cell lines.
%
% EXPCON         Whether to use expression constraints.  
%                Default is true.

%
envConstrain = 'Medium_MaxMinSign';
%

if ~exist('EXPCON', 'var')
    EXPCON = false;
end

protThreshDir = 'nci60prot_thresh';
micrThreshDir = 'nci60mRNA_thresh';

analysisLabel = '8DCLs';
CLs = {'MCF7', 'U251', 'COLO 205', 'CCRF-CEM', 'M14', 'NCI-H460', ...
       'SK-OV-3', 'PC-3'};

if ~exist('allCL', 'var')
    allCL = false;
end

if allCL
    analysisLabel = 'AllCLs';
    % !!! need to load CLs from NCI60_labels.csv
end

if ~exist('analysisOnly', 'var')
    analysisOnly = false;
end

if ~EXPCON
    analysisLabel = [analysisLabel '_noEXPCON'];
end

if exist('paramDir', 'var')
    if length(paramDir) > 0        
        analysisLabel = [analysisLabel '_' paramDir];
    end
end


% In case other constraints are added later
consString = '';
if exist('envConstrain', 'var')
    if length(envConstrain) > 0
        consString = [envConstrain '_'];
    end
end

% Run Simulations
if ~analysisOnly
    if ~exist('paramDir', 'var') || strcmp(paramDir, '')
        runMultiPerturbtion(model, protThreshDir, CLs, envConstrain, ...
                            analysisLabel, EXPCON);
        runMultiPerturbtion(model, micrThreshDir, CLs, envConstrain, ...
                            analysisLabel, EXPCON);
    else
        % make a parfor here and do a CL at a time
        parfor i = 1:length(CLs)
            runComparisonScript(model, 'FALCON', [protThreshDir '/' paramDir], ...
                                envConstrain, CLs(i), analysisLabel, EXPCON);
            runComparisonScript(model, 'FALCON', [micrThreshDir '/' paramDir], ...
                                envConstrain, CLs(i), analysisLabel, EXPCON);
        end
    end
end
falconProtOutDir = ['FALCON_' analysisLabel '_simresults_' consString protThreshDir '/'];
falconMicrOutDir = ['FALCON_' analysisLabel '_simresults_' consString micrThreshDir '/'];


% Do Analysis. colOrder = [3 1 2] (last column (3) specifies nan or 0).
analyzeMultiPerturbation(model, falconProtOutDir, [3 1 2]);
analyzeMultiPerturbation(model, falconMicrOutDir, [3 1 2]);

%[meanVals, stdVals] = analyzeMultiPerturbation_ErrorBars(...
%    analysisDir, suffix, splitLabels)