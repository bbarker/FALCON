function [v_solirrev v_solrev]=runFBA(model)
[modelIrrev,matchRev,rev2irrev,irrev2rev]=convertToIrreversible(model);

modelIrrev=changeObjective(modelIrrev,'biomass_reaction');
solutionStructIrrev=optimizeCbModel(modelIrrev,'max');
v_solirrev=solutionStructIrrev.x;

v_solrev=zeros(length(model.rxns),1);
for j=1:length(irrev2rev)
    irrevrxnname=modelIrrev.rxns{j};
    if~isempty(regexp(irrevrxnname,'_b$'))
        v_solrev(irrev2rev(j))=v_solrev(irrev2rev(j))-v_solirrev(j);
    elseif~isempty(regexp(irrevrxnname,'_f$'))
        v_solrev(irrev2rev(j))=v_solrev(irrev2rev(j))+v_solirrev(j);
    elseif~isempty(regexp(irrevrxnname,'_r$'))
        v_solrev(irrev2rev(j))=-v_solirrev(j);
    else
        v_solrev(irrev2rev(j))=v_solirrev(j);
    end
end
end