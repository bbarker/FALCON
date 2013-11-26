function modelNew=removeNonUsedGenes(model)
%Updating genes after reaction removal process finished
genes=unique(model.genes);
model.rules=[];
model.rxnGeneMat=[];
model.genes=[];
for m=1:numel(model.rxns)
    model = changeGeneAssociation(model,model.rxns{m,1},model.grRules{m,1},genes,genes);
end
modelNew=model;
