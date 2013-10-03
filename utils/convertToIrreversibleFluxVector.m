function v_solirrev=convertToIrreversibleFluxVector(modelIrrev,irrev2rev,v_solrev)
v_solirrev=zeros(length(modelIrrev.rxns),1);
for j=1:length(irrev2rev)
    irrevrxnname=modelIrrev.rxns{j};
    if (~isempty(regexp(irrevrxnname,'_b$')))
        if(v_solrev(irrev2rev(j))<0)
            v_solirrev(j)=-v_solrev(irrev2rev(j));
	end
    elseif (~isempty(regexp(irrevrxnname,'_f$')))
        if(v_solrev(irrev2rev(j))>0)
	    v_solirrev(j)=v_solrev(irrev2rev(j));
	end
    elseif ~isempty(regexp(irrevrxnname,'_r$'))
        v_solirrev(j)=-v_solrev(irrev2rev(j));
    else
        v_solirrev(j)=v_solrev(irrev2rev(j));
    end
end
end