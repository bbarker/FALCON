function [meanVals, stdVals] = analyzeMultiPerturbation_ErrorBars(...
    analysisDir, suffix, splitLabels)

%INPUT
% analysisDir
%
% suffix
%
%OPTIONAL INPUTS
% splitLabels    cell string of labels: assumes the first half of
%                each file represents one parameter, and the second
%                half the other parameter. 
% 

fileList = dir([analysisDir '/*' suffix '.csv']);
dataMat = [];

for i = 1:length(fileList)
    dataMat(i,:) = dlmread([analysisDir '/' fileList(i).name]);    
end

ncols = size(dataMat, 2);
meanVals = zeros(1, ncols);
stdVals = zeros(1, ncols);

%for i = 1:ncols
%    meanVals(i) = mean(dataMat(~isnan(dataMat(:, i)), i));
%    stdVals(i) = std(dataMat(~isnan(dataMat(:, i)), i));
%end

% Try to see if stability is the issue by checking
% several similar simulation thresholds
winRad = 2;
for i = (winRad + 1):(ncols-winRad)
    maxVals = zeros(length(fileList), 1);
    irng = i-winRad:i+winRad;
    for j = 1:length(fileList)
        maxVals(j) = max(dataMat(j,irng));
    end
    meanVals(i) = mean(maxVals(~isnan(maxVals)));
    stdVals(i) = std(maxVals(~isnan(maxVals)));
end

%if exist('splitLabels', 'var')
%end