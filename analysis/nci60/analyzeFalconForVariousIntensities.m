function dist = analyzeFalconForVariousIntensities (recMod, fileName, gene, rc, intenst, rxnNames, printPrefix, nV)
% This function looks at the change in flux for a reaction(s) when the 
% gene intensity value varies using Falcon. 

%INPUTS
% recMod is the reversible model (human recon 2 or yeast)
%
% fileName is the file to analyzed (ex. '786_O.csv')
%
% gene is a string of the gene of interest
%
% rc is regularization constant on fluxes
%
% intenst is a vector containing all gene intensities to be analyzed
%
% rxnNames is a cell array containing up to 3 reaction numbers 
%   ex. model.rxnNames{500}
%
% printPrefix is the name of the file where the graph will be printed (.jpg)
%
%OPTIONAL INPUT
% nV if it exists, is a rxnNames entry that will be used to normalize
% plotted fluxes.
%
%OUTPUTS
% dist is an array. 1st column is the gene intensity. Last column is 
% sum of all fluxes. Middle columns are fluxes of rxns of interest.
%
% Narayanan Sadagopan 10/12/13
% Brandon Barker      11/06/13  Added parallelization and minor changes.

rxnOfInt = cellfun(@(x) find(strcmp(recMod.rxnNames, x)), rxnNames);

EXPCON = false;
dist = zeros(length(intenst),length(rxnOfInt)+2);
ind = 0;

%create temporary file
fileID=fopen(fileName,'r');
c1 = textscan(fileID, '%s\t%s\t%s\n',1);
C = textscan(fileID, '%s\t%s\t%s\n');
cell2csv('tempFileForAnalyzeFalcon.csv', [c1{:}; C{:}], '\t', 2000); 
fclose(fileID);

%find gene of interest
for x = 1:length(C{1,1})
    if strcmp(gene, C{1,1}(x))
        ind = x;
        break;
    end
end
 
if (ind==0)
    disp('Gene not in model');
    return;
end

parfor x = 1:length(intenst)
    CC = C;
    CC{1,2}{ind} = num2str(intenst(x)); 
    cell2csv('tempFileForAnalyzeFalcon.csv', [c1{:}; CC{:}], '\t', 2000); 
    [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EXPCON, 0);
    dist(x, :) = [intenst(x) vRev(rxnOfInt)' norm(vRev, 1)];
    disp(dist(x, :));
end

delete('tempFileForAnalyzeFalcon.csv');

%plot
figure;
set(gca, 'FontSize', 36);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0.2 0.2 12 10]);
hold all
plot(dist(:, 1), dist(:, 2) ./ dist(:, end), 'b-o', 'MarkerSize', 10);
% why does this not work???? :
%hLeg = legend(rxnNames);
hLeg = legend(sprintf(rxnNames{1}));
if (length(rxnOfInt) == 2)
    plot(dist(:, 1), dist(:, 3)./ dist(:, end), 'g-*');
    hLeg = legend(sprintf(rxnNames{1}), sprintf(rxnNames{2}))
elseif (length(rxnOfInt) == 3)
    plot(dist(:, 1), dist(:, 3) ./ dist(:, end), 'g-*');
    plot(dist(:, 1), dist(:, 4) ./ dist(:, end), 'r-s');
    hLeg = legend(sprintf(rxnNames{1}), sprintf(rxnNames{2}), ...
        sprintf(rxnNames{3}));
end
set(hLeg, 'FontSize', 18);
title('Flux vs Intensity','FontSize', 36);
xlabel(sprintf('Intensity of Gene %s', gene), ...
    'FontSize', 32);
ylabel('Flux', 'FontSize', 32);
grid on;

print('-dpng', [printPrefix '.png'], '-r150');
print('-depsc2', [printPrefix '.eps'], '-r600');

set(hLeg, 'Location', 'Best');
print('-dpng', [printPrefix '_Best.png'], '-r150');
print('-depsc2', [printPrefix '_Best.eps'], '-r600');

set(hLeg, 'Location', 'BestOutside');
print('-dpng', [printPrefix '_BestOut.png'], '-r150');
print('-depsc2', [printPrefix '_BestOut.eps'], '-r600');

set(hLeg, 'Location', 'NorthEastOutside');
print('-dpng', [printPrefix '_NEOut.png'], '-r150');
print('-depsc2', [printPrefix '_NEOut.eps'], '-r600');

close(figure);