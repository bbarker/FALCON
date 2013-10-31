function pyrRatioPDH_PDC = customAnalyses(model, fluxFileDir)
%
% Do some specific analysis for each cell line in the 
% fluxFileDir

fluxFileList = dir([fluxFileDir '/*out']);
fluxFileList = struct2cell(fluxFileList);
% This assumes all the subdirectories have
% the same file names.
fluxFileList = sort(fluxFileList(1,:));
numCellLines = length(fluxFileList);

nrxns = length(model.rxns);
V_solrev = zeros(numCellLines, nrxns); 
for i = 1:numCellLines
    fluxFile = [fluxFileDir '/' fluxFileList{i}];
    [v_solex v_solrev v_solirrev] = readV_solFile(fluxFile);
    V_solrev(i, :) = v_solrev;
end

[cellLinesArray jainMetsArray coreTable ...
    FVAVminArray FVAVmaxArray] = readJainTable();

nci60filenames = cellfun(@convertExpressionFileName, cellLinesArray, ...
                         'UniformOutput', false);

%
% Calculate Warburg Effect and correlate with experimental lactate level.
%

% Get PDH and PDC indices in rev model.
PDH = find(strcmp(model.rxnNames, 'pyruvate dehydrogenase'));
LDH = find(strcmp(model.rxnNames, 'L-lactate dehydrogenase'))
%LDH = find(strcmp(model.rxnNames, ...
%          'L-Lactate dehydrogenase, cytosolic/mitochondrial'));

LACcore = find(strcmp(jainMetsArray, 'lac_D/lac_L'));

pyrRatioLDH_PDH = zeros(1, numCellLines);
for i = 1:numCellLines
    coreIdx = find(strcmp(nci60filenames, strrep(fluxFileList{i}, 'out', '')));
    if ~numel(coreIdx)
        fFile = fluxFileList{i}
        fFileBase = strrep(fluxFileList{i}, 'out', '')
    end
    disp([V_solrev(i, LDH)  V_solrev(i, PDH) coreTable(LACcore, coreIdx)])
    pyrRatioPDH_PDC(i) = sum(V_solrev(i, LDH)) / sum(V_solrev(i, PDH));
    %pyrRatioPDH_PDC
end

%nci60filenames