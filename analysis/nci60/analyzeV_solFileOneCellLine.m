function subsetsToStats = analyzeV_solFileOneCellLine(inputFI)
%
% Needs a spot of commenting in the code too
%
%

%%%% *** Settings ***
%
% Number of fields in subsetsToStats
nFields = 13;
%
%%%% ***    ***   *** 

[cellLinesArray jainMetsArray coreTable ...
    FVAVminArray FVAVmaxArray] = readJainTable();
subsetsToStats = containers.Map;
% !!! This shouldn't be the mean of the whole Table...
%coreTable = mean(coreTable);
[sortedCoreTableCol sortedCoreTableColIdxs] = sort(coreTable);
v_solExc = readV_solFile(inputFI);
sortedV_solExc = v_solExc(sortedCoreTableColIdxs);
%[sortedCoreTableCol sortedV_solExc]
sortedFVAVmaxArray = FVAVmaxArray(sortedCoreTableColIdxs);
sortedFVAVminArray = FVAVminArray(sortedCoreTableColIdxs);
for i = 1:length(sortedCoreTableColIdxs)
    statsArray = zeros(nFields,1);

    uptakeTruePos=0;
    uptakeFalseNeg=0;
    releaseTruePos=0;
    releaseFalseNeg=0;
    includedIdxs=[];
    for j=1:i
	if(sortedCoreTableCol(j) > 0 && sortedFVAVmaxArray(j) == 0)
	    statsArray(12) = statsArray(12) + 1;
	elseif(sortedCoreTableCol(j) < 0 && sortedFVAVmaxArray(j) == 0)
	    statsArray(12) = statsArray(12) + 1;
	elseif(sortedCoreTableCol(j) > 0)
	    statsArray(9) = statsArray(9) + 1;
	    includedIdxs(end + 1) = j;
	    if(sortedV_solExc(j) > 0)
		releaseTruePos = releaseTruePos + 1;
	    else
		releaseFalseNeg = releaseFalseNeg + 1;
	    end
	elseif(sortedCoreTableCol(j) < 0)
	    statsArray(8) = statsArray(8) + 1;
	    includedIdxs(end + 1) = j;
	    if(sortedV_solExc(j) < 0)
		uptakeTruePos = uptakeTruePos + 1;
	    else
		uptakeFalseNeg = uptakeFalseNeg + 1;
	    end
	else
	    statsArray(10) = statsArray(10) + 1;
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
	statsArray(1) = corr(Vest, Vexp, 'type', 'Pearson');
	statsArray(2) = corr(Vest, Vexp, 'type', 'Spearman');
	statsArray(3) = corr(Vest, Vexp, 'type', 'Kendall');
	statsArray(4) = Vest' * Vexp / (norm(Vest) * norm(Vexp));
	statsArray(5) = mean(Vest - Vexp);
	statsArray(6) = uptakeTruePos / (uptakeTruePos + uptakeFalseNeg);
	statsArray(7) = releaseTruePos / (releaseTruePos + releaseFalseNeg);
    end
    subsetsToStats(num2str(i))=statsArray;
end

lengthSubsetsToStats = length(keys(subsetsToStats));
subsetsToStats('Average') = zeros(13, 1);
subsetsToStats('Min') = subsetsToStats('1');
subsetsToStats('Max') = subsetsToStats('1');
numNotNaNArray = zeros(13, 1);
for j = 1:lengthSubsetsToStats
    subsetsToStats('Min') = [subsetsToStats('Min') ...
                             subsetsToStats(num2str(j))];
    subsetsToStats('Max') = [subsetsToStats('Max') ...
                             subsetsToStats(num2str(j))];

    jthStatsArray = subsetsToStats(num2str(j));
    newAverageStatsArray = subsetsToStats('Average');
    for k = 1:length(jthStatsArray)
        if ~isnan(jthStatsArray(k))
	  newAverageStatsArray(k) = newAverageStatsArray(k) + ...
                                    jthStatsArray(k);
	  numNotNaNArray(k) = numNotNaNArray(k) + 1;
        end
    end
    subsetsToStats('Average') = newAverageStatsArray;
end

subsetsToStats('Average') = subsetsToStats('Average') ./ numNotNaNArray;
subsetsToStats('Max') = max(subsetsToStats('Max'), [], 2);
subsetsToStats('Min') = min(subsetsToStats('Min'), [], 2);
