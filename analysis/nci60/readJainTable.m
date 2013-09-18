function [celllinesarray metsarray coretable FVAvminarray FVAvmaxarray] = ...
    readJainTable()
%INPUT
% Assumes Jain's table is in the running directory
%
%OUTPUT
%
%
warning('off', 'MATLAB:xlsread:ActiveX');
[excnumarray exctextarray raw] = xlsread( ... 
    ['Supp Table 3 A community-driven global reconstruction ' ... 
     'of human metabolism 95.xls']);
[height width] = size(excnumarray);
coretable1 = excnumarray(8:98, 8:width);
celllinesarray1 = exctextarray(9, 10:2:128);
coretable = [];
celllinesarray = {};
for i = 1:length(celllinesarray1)
    celllinesarray{end + 1} = celllinesarray1{i};
    coretable(:, end + 1) = mean(coretable1(:, (i*2-1):i*2), 2);
end
metsarray = exctextarray(10:100, 2);
FVAvminarray = excnumarray(8:98, 1);
FVAvmaxarray = excnumarray(8:98, 3);
end

