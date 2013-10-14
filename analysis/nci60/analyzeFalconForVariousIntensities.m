function dist = analyzeFalconForVariousIntensities (recMod, fileName, rc, intenst, rxnOfInt, printFile)
% This function looks at the change in flux for a reaction when the 
% gene intensity value varies using Falcon. For this function to work,
% the gene of interest who's intensity values are being changed MUST BE
% THE FIRST GENE in the file after the header

%INPUTS
%
% recMod is the reversible human recon 2 model
% fileName is the tissue file to analyzed (ex. '786_O.csv')
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


%
EXPCON = true;
%

dist = zeros(length(intenst),length(rxnOfInt)+1);
count = 1;

%create temporary file
fileID=fopen(fileName,'r');
c1 = textscan(fileID, '%s %s %s',1);
C = textscan(fileID, '%f %f %f');
hdr = {'gene', 'mean', 'var'}; 
txt=sprintf('%s\t', hdr{:});
txt(end) = '';
dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
dlmwrite('tempFileForAnalyzeFalcon.csv', C, '-append', 'delimiter', ...
         '\t', 'precision', 10);
fclose(fileID);

if (length(rxnOfInt)==1)
    %modify gene intensity and run Falcon
    for x = 1:length(intenst)

        %make change and rewrite file
        C{1,2}(1) = intenst(x); 
        dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
        dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter', ...
	     '\t','precision',10);

        %run Falcon
        [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EXPCON, 0);
        dist(count,1) = intenst(x);
        dist(count,2) = vRev(rxnOfInt(1));
        count = count + 1;
    end
elseif (length(rxnOfInt)==2)
    for x = 1:length(intenst)
        C{1,2}(1) = intenst(x); 
        dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
        dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter', ...
	     '\t','precision',10);
        [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EXPCON, 0);
        dist(count,1) = intenst(x);
        dist(count,2) = vRev(rxnOfInt(1));
        dist(count,3) = vRev(rxnOfInt(2));
        count = count + 1;
    end
else
    for x = 1:length(intenst)
        C{1,2}(1) = intenst(x); 
        dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
        dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter', ...
	     '\t','precision',10);
        [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc, EPCON, 0);
        dist(count,1) = intenst(x);
        dist(count,2) = vRev(rxnOfInt(1));
        dist(count,3) = vRev(rxnOfInt(2));
        dist(count,4) = vRev(rxnOfInt(3));
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
end
if (length(rxnOfInt)==3)
    plot(dist(:,1),dist(:,4),'r-o');
    legend(sprintf('reaction %d', rxnOfInt(1)),sprintf('reaction %d', ...
        rxnOfInt(2)),sprintf('reaction %d', rxnOfInt(3)),'Location', ...
        'NorthEastOutside');
end
title('Flux vs Intensity');
xlabel(sprintf('intensity of the gene %d', C{1,1}(1)));
ylabel(sprintf('flux for reaction %d', rxnOfInt));
grid on;
print('-dpng',printFile);
close(figure);