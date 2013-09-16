function runComparisonScript(model, method, expFileDir, envConstrain, CL)
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
% method:    one of: 
%                'FBA'
%                'LMOMA' (Linear MoMA)
%                'FALCON'
% 
% expFileDir    Relative path to top level directory containing
%               cell-line expression files. These files are
%               tab-delimited file with a head for the columns:
%               gene (entrez gene id), mean (expression value,
%               and standard deviation (of expression).
%
%
%OPTIONAL INPUTS
% envConstrain   'medium', 'core', or 'core_med':
%                whether or not to constraints that are 
%                medium-based CoRe-based, or both.
%
% CL             If nonempty or exists, should be the name
%                of a single cell-line to run.
%
%OUTPUT      a file in 'outputDir' (defined below) for each
%            cell line.
%
% Yiping Wang    08/??/13 
% Branon Barker  09/15/13   Changed output naming;
%                           added option for single cell line run.


mediumExcIdxs = loadMediumExcIdxs(model);
[celllinesarray jainMetsArray coretable] = readJainTable();
jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model);

%In case we only want to run a single cell line:
if exist('CL', 'var')
  clIdx = find(strcmp(celllinesarray, CL));
  celllinesarray = celllinesarray(clIdx);
end


for i = 1:length(celllinesarray)  
    consString = '';
    if exist('envConstrain', 'var')
        if length(envConstrain) > 0
            consString = [envConstrain '_'];
        end
    end
    modelToRun = model;
    if exist('envConstrain', 'var')
	if strcmp(envConstrain, 'medium')
	    modelToRun = constrainMediumExc(model, coretable(:,i));
	end
	if strcmp(envConstrain, 'core')
	    modelToRun = constrainCOREExc(model, coretable(:,i));
	end
	if strcmp(envConstrain, 'core_med')
	    modelToRun = constrainCOREExc(model, coretable(:,i));
	    modelToRun = constrainMediumExc(model, coretable(:,i));
	end
    end
    expressionFile = convertExpressionFileName(celllinesarray{i});
    outputDir = [method '_simresults_' consString expFileDir '/'];
    mkdir(outputDir);
    outputFile = [outputDir expressionFile 'out'];
    outputFI = fopen(outputFile, 'w');
    v_solirrev = [];
    v_solrev = [];
    if(strcmp(method, 'FALCON'))    
        expressionFileLoc = [expFileDir '/' expressionFile '.csv']
        [v_solirrev v_solrev] = runFalcon(modelToRun, ...
                                          expressionFileLoc, .1, .01);
    elseif(strcmp(method, 'LMOMA'))
        [v_solirrev v_solrev] = runLinearMOMA(model, ...
                                              rxnValues, rxnList);
    elseif(strcmp(method, 'FBA'))
        [v_solirrev v_solrev] = runFBA(model);
    else
        disp(['Error: unknown method ' method]);
	exit(-1);
    end

    fprintf(outputFI, 'All lower and upper bounds:\n');
    for j = 1:length(modelToRun.rxns)
        if(sum(mediumExcIdxs == j) ~= 0)
            fprintf(outputFI, '%s\t%20.15f\t%20.15f\n', ... 
                    modelToRun.rxns{j}, modelToRun.lb(j), ...
                    modelToRun.ub(j));
	else
	    for k = 1:length(jainMetsArray)
	        kthExcIdxs = jainMetsToExcIdxs(jainMetsArray{k});
		if(sum(kthExcIdxs == j) ~= 0)
		    fprintf(outputFI, '%s\t%20.15f\t%20.15f\n', ...
                            modelToRun.rxns{j}, ...
                            modelToRun.lb(j), modelToRun.ub(j));
		end
	    end
        end
    end

    fprintf(outputFI, 'All fluxes from v_solirrev:\n');
    for j = 1:length(v_solirrev)
        fprintf(outputFI, '%d\t%20.15f\n', j, v_solirrev(j));
    end
    
    fprintf(outputFI, 'All fluxes from v_solrev:\n');
    for j = 1:length(v_solrev)
        fprintf(outputFI, '%d\t%20.15f\n', j, v_solrev(j));
    end

    fprintf(outputFI, 'All fluxes from v_solex:\n');
    for j = 1:length(jainMetsArray)
	jthExcIdxs = jainMetsToExcIdxs(jainMetsArray{j});
        fprintf(outputFI, '%s\t%20.15f\n', jainMetsArray{j}, ...
                sum(v_solrev(jthExcIdxs)));
    end
    fclose(outputFI);
    end
end