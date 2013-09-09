function runComparisonScript(model,method,expFileDir)
mediumExcIdxs=loadMediumExcIdxs(model);
[celllinesarray jainMetsArray coretable]=readJainTable();
jainMetsToExcIdxs=loadJainMetsToExcIdxs(jainMetsArray,model);

for i=1:length(celllinesarray)  
    expressionFile=convertExpressionFileName(celllinesarray{i});
    outputDir = ['sims_' expFileDir '/'];
    mkdir(outputDir);
    outputFile=[outputDir expressionFile 'out'];
    outputFI=fopen(outputFile,'w');
    v_solirrev=[];
    v_solrev=[];
    %modelToRun=constrainCOREExc(model,coretable(:,i));
    %modelToRun=constrainMediumExc(model,coretable(:,i));
    modelToRun=model;    
    if(strcmp(method,'falcon'))    
        expressionFileLoc=[expFileDir '/' expressionFile '.csv']
        [v_solirrev v_solrev]=runFalcon(modelToRun,expressionFileLoc,.1,.01);
    elseif(strcmp(method,'linear MOMA'))
        [v_solirrev v_solrev]=runLinearMOMA(model,rxnValues,rxnList);
    elseif(strcmp(method,'FBA'))
        [v_solirrev v_solrev] = runFBA(model);
    else
        disp(['Error: unknown method ' method]);
	exit(-1);
    end

    fprintf(outputFI,'All lower and upper bounds:\n');
    for j=1:length(modelToRun.rxns)
        if(sum(mediumExcIdxs==j)~=0)
            fprintf(outputFI,'%s\t%20.15f\t%20.15f\n',modelToRun.rxns{j},modelToRun.lb(j),modelToRun.ub(j));
	else
	    for k=1:length(jainMetsArray)
	        kthExcIdxs=jainMetsToExcIdxs(jainMetsArray{k});
		if(sum(kthExcIdxs==j)~=0)
		    fprintf(outputFI,'%s\t%20.15f\t%20.15f\n',modelToRun.rxns{j},modelToRun.lb(j),modelToRun.ub(j));
		end
	    end
        end
    end

    fprintf(outputFI,'All fluxes from v_solirrev:\n');
    for j=1:length(v_solirrev)
        fprintf(outputFI,'%d\t%20.15f\n',j,v_solirrev(j));
    end
    
    fprintf(outputFI,'All fluxes from v_solrev:\n');
    for j=1:length(v_solrev)
        fprintf(outputFI,'%d\t%20.15f\n',j,v_solrev(j));
    end

    fprintf(outputFI,'All fluxes from v_solex:\n');
    for j=1:length(jainMetsArray)
	jthExcIdxs=jainMetsToExcIdxs(jainMetsArray{j});
        disp(j);
        disp(jainMetsArray{j}); 
        disp(jthExcIdxs);
        disp(length(v_solrev));
        disp('===================')
        fprintf(outputFI,'%s\t%20.15f\n',jainMetsArray{j},sum(v_solrev(jthExcIdxs)));
    end
    fclose(outputFI);
    end
end