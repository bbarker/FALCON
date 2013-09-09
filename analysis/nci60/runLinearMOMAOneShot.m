function [v_solirrev v_solrev]=runLinearMOMAOneShot(model,rxnValues,rxnList)
[modelIrrev,matchRev,rev2irrev,irrev2rev]=convertToIrreversible(model);
WTflux=zeros(length(model.rxns),1);
for i=1:length(rxnList)
    WTflux(rxnList(i))=rxnValues(i);
end
WTflux=convertToIrreversibleFluxVector(modelIrrev,irrev2rev,WTflux)
[solutionDel,solutionWT]=linearMOMAOneShot(modelIrrev,WTflux,0,0,1,zeros(length(modelIrrev.rxns),1),rxnList);
v_solirrev=solutionDel.x
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