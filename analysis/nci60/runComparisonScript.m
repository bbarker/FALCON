function runComparisonScript(model, method, expFileDir, envConstrain, ...
    CLs, addLabel, EXPCON)
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
% CLs            If nonempty or exists, should be the name
%                of list the cell-lines to run, as found in
%                the CoRe table.
%
% addLabel       Label for any changes made to the model 
%                to be used in naming the output directory.
%
% EXPCON         Whether to use expression constraints.  
%                Default is true.
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

sz_cla = size(celllinesarray)
%In case we only want to run a single cell line:
if exist('CLs', 'var')
    clIdx = [];
    for i = 1:length(CLs)
        clIdx = [clIdx find(strcmp(celllinesarray, CLs{i}))]
    end
end
celllinesarray = celllinesarray(clIdx);
sz_cla = size(celllinesarray)

if ~exist('EXPCON', 'var')
    EXPCON = false;
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
    modLabel = '';
    if exist('addLabel', 'var')
        if length(addLabel) > 0
            modLabel = ['_' addLabel];
        end
    end
    outputDir = [method modLabel '_simresults_' consString expFileDir '/'];
    mkdir(outputDir);
    outputFile = [outputDir expressionFile 'out'];
    outputFI = fopen(outputFile, 'w');
    v_solirrev = [];
    v_solrev = [];
    if(strcmp(method, 'FALCON'))    
        expressionFileLoc = [expFileDir '/' expressionFile '.csv']
        [v_solirrev v_solrev] = runFalcon(modelToRun, ...
                                          expressionFileLoc, .01, EXPCON);
    elseif(strcmp(method, 'LMOMA'))
        rxnValues = [];
        rxnList = [];
	for j = 1:length(jainMetsArray)
	    jthExcIdxs = jainMetsToExcIdxs(jainMetsArray{j});
            for k = 1:length(jthExcIdxs)
                rxnList(end + 1) = jthExcIdxs(k);
                % Assume the maximum value (usually glucose) is always attainable
                % at unity:
                rxnValues(end + 1) = coretable(j, i) / max(coretable(:, i));
                % Everything appears to check out:
                % disp([celllinesarray{i} '    ' jainMetsArray{j}]);
		% disp(jthExcIdxs(k));
		% disp(coretable(j, i));
            end
	end
        [v_solirrev v_solrev] = runLinearMOMAOneShot(model, ...
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