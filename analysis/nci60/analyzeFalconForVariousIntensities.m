function dist = analyzeFalconForVariousIntensities (recMod, fileName, gene, rc, intenst, rxnNames, printFile)
% This function looks at the change in flux for a reaction(s) when the 
% gene intensity value varies using Falcon. 

%INPUTS
%
% recMod is the reversible human recon 2 model
% fileName is the tissue file to analyzed (ex. '786_O.csv')
% gene is the number of the gene of interest
% rc is regularization constant on fluxes
% intenst is a vector containing all gene intensities to be analyzed
% rxnOfInt is a vector containing up to 3 reaction numbers
% printFile is the name of the file where the graph will be printed (.jpg)
%
%OUTPUTS
% dist is an array containing 2 columns. 1st column is the gene 
% intensity and the 2nd column is the reaction flux. 
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
hdr = {'gene', 'mean', 'var'}; 
txt=sprintf('%s\t', hdr{:});
txt(end) = '';
sz_txt = size(c1)
sz_C = size(C)
sz_C = size( C{1,2})
sz_C = size( C{1, 3})
%Cadj = {C{1}; C{2}; C{3}};
%sz_C = size(Cadj)

cell2csv('tempFileForAnalyzeFalcon.csv', [c1{:}; C{:}], '\t', 2000); 
%dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
%dlmwrite('tempFileForAnalyzeFalcon.csv', C, '-append', 'delimiter', ...
%         '\t', 'precision', 10);
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
        %dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
        %dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter', ...
        %     '\t','precision',10);

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
        %dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
        %dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter', ...
        %     '\t','precision',10);
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
        %dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
        %dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter', ...
        %     '\t','precision',10);
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
plot(dist(:,1),dist(:,2),'b-o');
legend(sprintf('reaction %d', rxnOfInt(1)),'Location','NorthEastOutside');
if (length(rxnOfInt)==2)
    plot(dist(:,1),dist(:,3),'g-o');
    legend(sprintf('reaction %d', rxnOfInt(1)),sprintf('reaction %d', ...
        rxnOfInt(2)),'Location','NorthEastOutside');
elseif (length(rxnOfInt)==3)
    plot(dist(:,1),dist(:,3),'g-o');
    plot(dist(:,1),dist(:,4),'r-o');
    legend(sprintf('reaction %d', rxnOfInt(1)),sprintf('reaction %d', ...
        rxnOfInt(2)),sprintf('reaction %d', rxnOfInt(3)),'Location', ...
        'NorthEastOutside');
end
title('Flux vs Intensity');
xlabel(sprintf('intensity of the gene with Entrez ID %d', ind));
ylabel('flux');
grid on;
print('-dpng',printFile);
close(figure);