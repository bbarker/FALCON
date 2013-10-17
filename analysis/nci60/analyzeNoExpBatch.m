function analyzeNoExpBatch(expFileDir)
%
% Used to analyze and summarize stats across all cell lines for 
% simulation methods that do not use expression, such as FBA and 
% Linear MoMA.


% How many ranked exo metabolites to consider
metDepth = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 20 30 40 50 60 70 80 90]; 

% Which stats to get and their labels.
statLabel = {'Pearson', 'Spearman', 'Kendall', 'CosineSim', 'L1Dist', ...
             'Sensitivity', 'UptakeSensitivity', 'ReleaseSensitivity'};

tab = sprintf('\t');

dirOut = [rstrtrim(expFileDir, '/') '_stats'];
mkdir(dirOut);
fluxFileList = dir([expFileDir '/*out']);
fluxFileList = struct2cell(fluxFileList);
fluxFileList = sort(fluxFileList(1,:));
numCellLines = length(fluxFileList);
subsetsToStatsCell = cell([1, numCellLines]);
parfor i = 1:numCellLines
    fluxFile = [expFileDir '/' fluxFileList{i}];
    subsetsToStatsCell{i} = analyzeV_solFileOneCellLine(fluxFile);
end

summaryFileName = [dirOut '/summary.csv'];
summaryFI = fopen(summaryFileName, 'w');
fprintf(summaryFI, '%s\t%s\t%s\n', 'stat', 'mean', 'std');

for k = 1:length(metDepth)
    dstr = num2str(metDepth(k));
    for s = 1:length(statLabel)
	skey = statLabel{s};
	clString = '';
	valString = '';
        valVec = [];
	for i = 1:numCellLines
            statsCell = subsetsToStatsCell{i}(dstr);
            clString = sprintf('%s\t%s', clString, fluxFileList{i});
	    valString = sprintf('%s\t%s', valString, ...
				num2str(statsCell(skey)));
	    valString = lstrtrim(valString, tab);
            valVec(end + 1) = statsCell(skey); 
        end
        fprintf(summaryFI, '%s\t%s\t%s\n', [statLabel{s} dstr], ...
            num2str(mean(valVec(~isnan(valVec)))), ...
            num2str(std(valVec(~isnan(valVec)))));
	fileName = [dirOut '/' statLabel{s} dstr '.csv'];
	flatAnalysisFI = fopen(fileName, 'w');
	fprintf(flatAnalysisFI, '%s\n', clString);
	fprintf(flatAnalysisFI, '%s\n', valString);
	fclose(flatAnalysisFI);
    end
end

fclose(summaryFI);