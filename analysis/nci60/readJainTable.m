function [celllinesarray metsarray coretable FVAvminarray FVAvmaxarray]=readJainTable()
[excnumarray exctextarray raw]=xlsread('Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
coretable1=excnumarray(8:98,8:width);
celllinesarray1=exctextarray(9,10:2:128);
coretable=[];
celllinesarray={};
for i=1:length(celllinesarray1)
    if(~strcmp(celllinesarray1{i},'MDA-MB-468') && ~strcmp(celllinesarray1{i},'RXF 393'))
        celllinesarray{end+1}=celllinesarray1{i};
	coretable(:,end+1)=coretable1(:,i*2-1);
	%coretable(:,end+1)=coretable1(:,i*2);
    end
end
metsarray=exctextarray(10:100,2);
FVAvminarray=excnumarray(8:98,1);
FVAvmaxarray=excnumarray(8:98,3);
end

