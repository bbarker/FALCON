function models = constrainCOREExcequential(modelOriginal, coreValues)
%This function takes an original model, usually rec2, and outputs a cell array of models
%with the bounds on certain reactions sequentially constrained, depending on the measured CORE values from Jain et al.

[cellLinesArray jainMetsArray coreTable]=readJainTable();
jainMetsToExcIdxs=loadJainMetsToExcIdxs(jainMetsArray,model);
mediumExcIdxs=loadMediumExcIdxs(model);
%sort CORE values with greatest absolute value first.
[junk coreValuesAbsSortInds]=sort(abs(coreValues),'descend');

%Constrained all exchange reactions not in the medium to have no uptake, lb=0.
modelConstrainedRelease=modelOriginal;
for i=1:length(modelConstrainedRelease.rxns)
    %if (~isempty(regexp(modelConstrainedRelease.rxns{i},'^EX(.)*\(e\)$')))
        %foundInJainMetsToExcIdxs=0;
	%for k=1:length(jainMetsArray)
	%    kthExcIdxs=jainMetsToExcIdxs(jainMetsArray{k});
	%    if(sum(kthExcIdxs==j)~=0)
	%	foundInJainMetsToExcIdxs=1;    
	%    end
	%end
    if (sum(mediumExcIdxs==i)==0)
        modelConstrainedRelease.lb(i)=0;
    end
    %end
end

%Sequentially set lbs of CORE exc rxns by narayan's scaling,
%starting with largest absolute value of CORE.
models{1}=modelConstrainedRelease;
for i=1:length(coreValuesAbsSortIdxs)
    modelConstrainedSeq=models{end};
    excIdxs=jainMetsToExcIdxs(jainMetsArray{coreValuesAbsSortIdxs(i)});
    modelChanged=1;
    for j=1:length(excIdxs)
        foundInJainMetsToExcIdxs=0;
	for k=1:length(jainMetsArray)
	    kthExcIdxs=jainMetsToExcIdxs(jainMetsArray{k});
	    if(sum(kthExcIdxs==j)~=0)
		foundInJainMetsToExcIdxs=1;    
	    end
	end
        if (foundInJainMetsToExcIdxs && coreValues(coreValuesAbsSortIdxs(i))<=0)
	    modelConstrainedSeq.lb(excIdxs(j))=coreValues(coreValuesAbsSortIdxs(i))/300;
	elseif (foundInJainMetsToExcIdxs && coreValues(coreValuesAbsSortIdxs(i))>=0)
            modelConstrainedSeq.ub(excrxninds(j))=coreValues(coreValuesAbsSortIdxs(i))/300;
	else
	    modelChanged=0;
        end
    end
    if(modelChanged)
        models{end+1}=modelConstrainedSeq;
    end
end