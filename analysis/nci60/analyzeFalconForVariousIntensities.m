function dist = analyzeFalconForVariousIntensities (recMod, fileName, gene, rc, intenst, rxnNames, printFile)
% This function looks at the change in flux for a reaction(s) when the 
% gene intensity value varies using Falcon. 

%INPUTS
%
% recMod is the reversible model (human recon 2 or yeast)
% fileName is the file to analyzed (ex. '786_O.csv')
% gene is a string of the gene of interest
% rc is regularization constant on fluxes
% intenst is a vector containing all gene intensities to be analyzed
% rxnNames is a cell array containing up to 3 reaction numbers 
%   ex. model.rxnNames{500}
% printFile is the name of the file where the graph will be printed (.jpg)
%
%OUTPUTS
% dist is an array. 1st column is the gene intensity. Last column is 
% sum of all fluxes. Middle columns are fluxes of rxns of interest.
%
% Narayanan Sadagopan 10/12/13

rxnOfInt = cellfun(@(x) find(strcmp(recMod.rxnNames, x)), rxnNames);

EXPCON = false;
dist = zeros(length(intenst),length(rxnOfInt)+2);
count = 1;
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

if (length(rxnOfInt)==1)
    %modify gene intensity and run Falcon
    for x = 1:length(intenst)

        %make change and rewrite file
        C{1,2}{ind} = num2str(intenst(x));
        cell2csv('tempFileForAnalyzeFalcon.csv', [c1{:}; C{:}], '\t', 2000); 

        %run Falcon
        [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EXPCON, 0);
        dist(count,1) = intenst(x);
        dist(count,2) = vRev(rxnOfInt(1));
        dist(count,3) = norm(vRev, 1);
        disp(dist(count,:));
        count = count + 1;
    end
elseif (length(rxnOfInt)==2)
    for x = 1:length(intenst)
        C{1,2}{ind} = num2str(intenst(x)); 
        cell2csv('tempFileForAnalyzeFalcon.csv', [c1{:}; C{:}], '\t', 2000); 
        [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EXPCON, 0);
        dist(count,1) = intenst(x);
        dist(count,2) = vRev(rxnOfInt(1));
        dist(count,3) = vRev(rxnOfInt(2));
        dist(count,4) = norm(vRev, 1);
        disp(dist(count,:));
        count = count + 1;
    end
else
    for x = 1:length(intenst)
        C{1,2}{ind} = num2str(intenst(x)); 
        cell2csv('tempFileForAnalyzeFalcon.csv', [c1{:}; C{:}], '\t', 2000); 
        [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EXPCON, 0);
        dist(count,1) = intenst(x);
        dist(count,2) = vRev(rxnOfInt(1));
        dist(count,3) = vRev(rxnOfInt(2));
        dist(count,4) = vRev(rxnOfInt(3));
        dist(count,5) = norm(vRev, 1);
        disp(dist(count,:));
        count = count + 1;
    end
end

delete('tempFileForAnalyzeFalcon.csv');

%plot
figure;
hold all
plot(dist(:,1),dist(:,2),'b-o','MarkerSize',10);
legend(sprintf(rxnNames{1}),'Location','NorthEastOutside');
if (length(rxnOfInt)==2)
    plot(dist(:,1),dist(:,3),'g-*');
    legend(sprintf(rxnNames{1}),sprintf(rxnNames{2}), ...
        'Location','NorthEastOutside');
elseif (length(rxnOfInt)==3)
    plot(dist(:,1),dist(:,3),'g-*');
    plot(dist(:,1),dist(:,4),'r-s');
    legend(sprintf(rxnNames{1}),sprintf(rxnNames{2}), ...
        sprintf(rxnNames{3}));
end
title('Flux vs Intensity','FontSize',20);
xlabel(sprintf('Intensity of the Gene with ID %s', gene), ...
    'FontSize',16);
ylabel('Flux','FontSize',16);
grid on;
print('-dpng',printFile);
close(figure);