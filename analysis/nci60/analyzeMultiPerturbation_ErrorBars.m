function [meanVals, stdVals] = analyzeMultiPerturbation_ErrorBars(...
    analysisDir, suffix, splitLabels)

%INPUT
% analysisDir    Output from analyzeMultiPerturbation.m
%
% suffix         Which analysis (e.g. 'Pearson5') to plot.
%
%OPTIONAL INPUTS
% splitLabels    cell string of labels: assumes the first half of
%                each file represents one parameter, and the second
%                half the other parameter. 
% 

figExt = '.png';


fileList = dir([analysisDir '/*' suffix '.csv']);
dataMat = [];

for i = 1:length(fileList)
    dataMat(i, :) = dlmread([analysisDir '/' fileList(i).name]);    
end

ncols = size(dataMat, 2);
meanVals = zeros(1, ncols);
stdVals = zeros(1, ncols);


if exist('splitLabels', 'var')
    1+1

else
     [meanVals, stdVals] = performPlotting(dataMat, ...
                               [analysisDir '/' suffix figExt]);
end


end % of analyzeMultiPerturbation_ErrorBars

function [mVals, sVals] = performPlotting(dataSubMat, outFileName)
% This function can be called just once, or one time for each
% splitLabels
%

nrows = size(dataSubMat, 1);
ncols = size(dataSubMat, 2);
mVals = zeros(1, ncols);
sVals = zeros(1, ncols);

%for i = 1:ncols
%    mVals(i) = mean(dataSubMat(~isnan(dataSubMat(:, i)), i));
%    sVals(i) = std(dataSubMat(~isnan(dataSubMat(:, i)), i));
%end

% Try to see if stability is the issue by checking
% several similar simulation thresholds

winRad = 2;
for i = (winRad + 1):(ncols - winRad)
    winMeans = zeros(nrows, 1);
    irng = (i - winRad):(i + winRad);
    for j = 1:nrows
        winMeans(j) = mean(dataSubMat(j, irng));
    end
    mVals(i) = mean(winMeans(~isnan(winMeans)));
    sVals(i) = std(winMeans(~isnan(winMeans)));
end

end % end of performPlotting