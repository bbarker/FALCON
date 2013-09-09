function models = constrainCOREExcequential(modeloriginal, coreValues)
%This function takes an original model, usually rec2, and outputs a cell array of models
%with the bounds on certain reactions sequentially constrained, depending on .

[celllinesarray jainMetsArray coretable]=readJainTable();
jainMetsToExcIdxs=loadJainMetsToExcIdxs(jainMetsArray,model);
mediumExcIdxs=loadMediumExcIdxs(model);
[junk abssortInds]=sort(abs(coreValues),'descend');

modelConstrainedCOREExc=modeloriginal;
for i=1:length(modelConstrainedCOREExc.rxns)
    if (~isempty(regexp(modelConstrainedCOREExc.rxns{i},'^EX(.)*\(e\)$')))
        foundInJainMetsToExcIdxs=0;
	for k=1:length(jainMetsArray)
	    kthExcIdxs=jainMetsToExcIdxs(jainMetsArray{k});
	    if(sum(kthExcIdxs==j)~=0)
		foundInJainMetsToExcIdxs=1;    
	    end
	end
        if (sum(mediumExcIdxs==i)==0 && ~foundInJainMetsToExcIdxs)
            modelConstrainedCOREExc.lb(i)=0;
        end
    end
end

models{1}=modelConstrainedCoreExc;
for i=1:length(abssortInds)
    modelsequnconstrained=models{end};
    excIdxs=jainMetsToExcIdxs(jainMetsArray{abssortInds(i)});
    modelChanged=0;
    for j=1:length(excIdxs)
        foundInJainMetsToExcIdxs=0;
	for k=1:length(jainMetsArray)
	    kthExcIdxs=jainMetsToExcIdxs(jainMetsArray{k});
	    if(sum(kthExcIdxs==j)~=0)
		foundInJainMetsToExcIdxs=1;    
	    end
	end
        if (foundInJainMetsToExcIdxs && coreValues(abssortInds(i))<=0)
	  modelsequnconstrained.lb(excIdxs(j))=coreValues(abssortInds(i))/300;
	  modelChanged=1;
	elseif (foundInJainMetsToExcIdxs && coreValues(abssortInds(i))>=0)
          modelsequnconstrained.ub(excrxninds(j))=coreValues(abssortInds(i))/300;
	  modelChanged=1;
	elseif (sum(mediumExcIdxs==i)==0 && coreValues(abssortInds(i))<=0)
	    modelsequnconstrained.lb(excrxninds(j))=coreValues(abssortInds(i))/300;
	    modelChanged=1;
        elseif (sum(mediumExcIdxs==i)==0 && coreValues(abssortInds(i))>=0)
	    modelsequnconstrained.ub(excrxninds(j))=coreValues(abssortInds(i))/300;
	    modelChanged=1;
        end
    end
    if(modelChanged)
        models{end+1}=modelsequnconstrained;
    end
end