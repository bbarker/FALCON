function ind2EX=findExRxns(model)
%custom implementaion of identifying exchange reactions in Recon1
ind2EX=~cellfun(@isempty,strfind(model.rxns,'EX_'));