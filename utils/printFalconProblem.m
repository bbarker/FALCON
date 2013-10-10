function printFalconProblem(rowLabels, colLabels, cnt, A, b, lb, ub, f, ...
                            csense, y)
%
% Used for debugging falcon.m. Normally not called.
% Could be extended to print additional problem info.
%

[nrows, ncols] = size(A);
nOutCols = ncols + 1 + 2;           %+2 for each additional column,
nOutRows = nrows + 8;               %so as to include space
outCell = cell(nOutRows, nOutCols);
outCell{1, 1} = '';

for i = 1:nrows
    for j = 1:ncols
        if A(i,j) ~= 0
            outCell{i+1, j+1} = num2str(A(i, j));
        else
            outCell{i+1, j+1} = '';
        end
    end
end

%disp(size(outCell))
%disp(size(A))
%disp(size(lb))
%disp(size(colLabels))
%disp(ncols)
%disp(nrows)
%disp(colLabels)

for i = 1:ncols
    outCell{1, i+1} = colLabels{i};
end
outCell{1, ncols + 1 + 1} = 'csens';
outCell{1, ncols + 1 + 2} = 'b';
for i = 1:nrows
    outCell{1 + i, ncols + 2} = csensePrint(csense(i));
    outCell{i + 1, ncols + 3} = num2str(b(i));
end

outCell{nrows + 1 + 2, 1} = 'LB';
outCell{nrows + 1 + 4, 1} = 'UB';
outCell{nrows + 1 + 6, 1} = 'Obj';
outCell{nrows + 1 + 8, 1} = 'y';
for i = 1:ncols
    outCell{nrows + 2, i + 1} = '';
    outCell{nrows + 3, i + 1} = num2str(lb(i));
    outCell{nrows + 4, i + 1} = '';
    outCell{nrows + 5, i + 1} = num2str(ub(i));
    outCell{nrows + 6, i + 1} = '';
    outCell{nrows + 7, i + 1} = num2str(f(i));
    outCell{nrows + 8, i + 1} = '';
    if numel(y) <= 1
        outCell{nrows + 9, i + 1} = '';
    else
        outCell{nrows + 9, i + 1} = num2str(y(i));
    end
end


for i = 1:nrows
    outCell{i+1, 1} = rowLabels{i};
end

cell2csv(['FalconProblem' num2str(cnt) '.csv'], outCell, ',', 2000);
end % of printFalconProblem

function reln = csensePrint(achar)
reln = 'ERROR';
if achar == 'L'
    reln = '<=';
elseif achar == 'E'
    reln = '=';
elseif achar == 'G'
    reln = '>=';
end
end % of csensePrint

