function subsetsToStats = analyzeV_solFileOneCellLine(inputFI, fluxCol)
%
% Needs a spot of commenting in the code too
%
%

%OPTIONAL INPUTS
% fluxCol    (1-60) - specifies an alternative CORE flux to compare to,
%                     primarily for the purpose of statistical analysis
%

%%%% *** Settings ***
%
% Number of fields in subsetsToStats
nFields = 13;
%
%%%% ***    ***   *** 

[cellLinesArray jainMetsArray coreTable ...
    FVAVminArray FVAVmaxArray] = readJainTable();

% Need to find index of file in coretable
inputFiSplit = strsplit(inputFI, '/');
fileName = inputFiSplit{end};
CLidx = -1;
if exist('fluxCol', 'var')
    CLidx = fluxCol;
else
    for i = 1:length(cellLinesArray)
	CL = convertExpressionFileName(cellLinesArray{i});
	empty = isempty(regexp(fileName, ['^' CL]));
	if ~empty
	    CLidx = i;
	    break;
	end
    end
end
if CLidx <= 0
    disp('Error: cell line not found.');
    return;
end



subsetsToStats = containers.Map;
[sortedCoreTableCol sortedCoreTableColIdxs] = sort(abs(coreTable(:, CLidx)), 1, 'descend');
sortedCoreTableCol = coreTable(sortedCoreTableColIdxs, CLidx);
v_solExc = readV_solFile(inputFI);
sortedV_solExc = v_solExc(sortedCoreTableColIdxs);
sortedFVAVmaxArray = FVAVmaxArray(sortedCoreTableColIdxs);
sortedFVAVminArray = FVAVminArray(sortedCoreTableColIdxs);
for i = 1:length(sortedCoreTableColIdxs)
    %statsCell = zeros(nFields,1);
    statsCell = containers.Map;
    statsCell('8') = 0;
    statsCell('9') = 0;
    statsCell('10') = 0;
    statsCell('12') = 0;

    uptakeTruePos=0;
    uptakeFalseNeg=0;
    releaseTruePos=0;
    releaseFalseNeg=0;
    includedIdxs=[];
    for j=1:i
	if(sortedCoreTableCol(j) > 0 && sortedFVAVmaxArray(j) == 0)
	    statsCell('12') = statsCell('12') + 1;
	elseif(sortedCoreTableCol(j) < 0 && sortedFVAVmaxArray(j) == 0)
	    statsCell('12') = statsCell('12') + 1;
	elseif(sortedCoreTableCol(j) > 0)
	    statsCell('9') = statsCell('9') + 1;
	    includedIdxs(end + 1) = j;
	    if(sortedV_solExc(j) > 0)
		releaseTruePos = releaseTruePos + 1;
	    else
		releaseFalseNeg = releaseFalseNeg + 1;
	    end
	elseif(sortedCoreTableCol(j) < 0)
	    statsCell('8') = statsCell('8') + 1;
	    includedIdxs(end + 1) = j;
	    if(sortedV_solExc(j) < 0)
		uptakeTruePos = uptakeTruePos + 1;
	    else
		uptakeFalseNeg = uptakeFalseNeg + 1;
	    end
	else
	    statsCell('10') = statsCell('10') + 1;
	    includedIdxs(end + 1) = j;
	end
    end
    uptakeFalsePos = releaseFalseNeg;
    uptakeTrueNeg = releaseTruePos;
    releaseFalsePos = uptakeFalseNeg;
    releaseTrueNeg = uptakeTruePos;

    if(~isempty(includedIdxs))
        %figure('Position', [300, 300, 1200, 500], 'Visible', 'off');
	%[trimmedcoretablecol trimmedv_solex] = ...
        %    trimOutliers(coretable(includedinds,i),v_solex(includedinds));
	%a=scatter(trimmedcoretablecol, trimmedv_solex);
	%set(gca,'Title',text('String',[expressionFile 'scatter']));
	%saveas(a,[inputfolder expressionFile 'scatter.png']);
	%statsarray(1)=corr(trimmedv_solex,trimmedcoretablecol,'type','Pearson');
        %statsarray(2)=corr(trimmedv_solex,trimmedcoretablecol,'type','Spearman');
        %statsarray(3)=corr(trimmedv_solex,trimmedcoretablecol,'type','Kendall');
        Vest = columnVector(sortedV_solExc(includedIdxs));
        Vexp = columnVector(sortedCoreTableCol(includedIdxs));
	statsCell('Pearson') = corr(Vest, Vexp, 'type', 'Pearson');
	statsCell('Spearman') = corr(Vest, Vexp, 'type', 'Spearman');
	statsCell('Kendall') = corr(Vest, Vexp, 'type', 'Kendall');
	statsCell('CosineSim') = Vest' * Vexp / (norm(Vest) * norm(Vexp));
	% Changing this to norm to account for cancellation.
	% statsCell(5) = mean(Vest - Vexp);
	statsCell('L1Dist') = norm(Vest - Vexp, 1);
        statsCell('Sensitivity') = (uptakeTruePos + releaseTruePos) / ...
                                   (uptakeTruePos + uptakeFalseNeg + ...
                                   releaseTruePos + releaseFalseNeg);
	statsCell('UptakeSensitivity') = uptakeTruePos / ...
                                          (uptakeTruePos + uptakeFalseNeg);
	statsCell('ReleaseSensitivity') = releaseTruePos / ... 
                                          (releaseTruePos + releaseFalseNeg);
    end
    subsetsToStats(num2str(i)) = statsCell;
end

statKeys = keys(subsetsToStats('1'));
lengthSubsetsToStats = length(keys(subsetsToStats));

subsetsToStats('Average') = subsetsToStats('1');
subsetsToStats('Min') = subsetsToStats('1');
subsetsToStats('Max') = subsetsToStats('1');
nnanCount = 0;
for i =  1:length(statKeys)
    avgCell = subsetsToStats('Average');
    minCell = subsetsToStats('Min');
    maxCell = subsetsToStats('Max');
    mysum = 0;
    mymax = -Inf;
    mymin = Inf;
    statsCell = subsetsToStats(num2str(i));
    mykey = statKeys{i};
    for j = 1:lengthSubsetsToStats
        statVal = statsCell(mykey);
        if ~isnan(statVal)
            mysum = mysum + statsCell(mykey);
            mymin = min(mymin, statsCell(mykey));
            mymax = max(mymax, statsCell(mykey));
        else
            nnanCount = nnanCount + 1;
        end
    end
    avgCell(mykey) = mysum / (lengthSubsetsToStats - nnanCount);
    minCell(mykey) = mymin;
    maxCell(mykey) = mymax;
    subsetsToStats('Average') = avgCell;
    subsetsToStats('Min') = minCell;
    subsetsToStats('Max') = maxCell;
end

