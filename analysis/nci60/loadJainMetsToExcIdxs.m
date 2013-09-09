function jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray,model)
jainMetsToExcIdxs=containers.Map;
for i=1:length(jainMetsArray)
    if(strcmp(jainMetsArray{i},'34hpp'))
        excrxnname='EX_34hpp';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    elseif(strcmp(jainMetsArray{i},'glc_D'))
        excrxnname='EX_glc(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    elseif(strcmp(jainMetsArray{i},'udpgal/udpg'))
        excrxnname1='EX_udpgal(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_udpg(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1(1) excrxnind2(1)];
    elseif(strcmp(jainMetsArray{i},'lac_D/lac_L'))
        excrxnname1='EX_lac_D(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_lac_L(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1(1) excrxnind2(1)];
    elseif(strcmp(jainMetsArray{i},'sbt_D'))
        excrxnname='EX_sbt-d(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    elseif(strcmp(jainMetsArray{i},'tyr_l'))
        excrxnname='EX_tyr_L(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    elseif(strcmp(jainMetsArray{i},'cit/icit'))
        excrxnname='EX_cit(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    elseif(strcmp(jainMetsArray{i},'N/A'))
        excrxnname='';
        excrxnind=0;
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    else
        excrxnname=strcat('EX_',strcat(jainMetsArray{i},'(e)'));
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind];
    end
end
end

