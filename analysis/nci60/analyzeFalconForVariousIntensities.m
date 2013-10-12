%This function looks at the change in flux for a reaction when the 
%gene intensity value varies using Falcon. For this function to work,
%the gene of interest who's intensity values are being changed MUST BE
%THE FIRST GENE in the file after the header
function dist = analyzeFalconForVariousIntensities (recMod, fileName, rc, intenst, rxnOfInt, printFile)

%INPUTS
    %recMod is the reversible human recon 2 model
    %fileName is the tissue file to analyzed (ex. '786_O.csv')
    %rc is regularization constant on fluxes
    %intenst is a vector containing all gene intensities to be analyzed
    %rxnOfInt is the number of the rxn where fluxes will be measured
    %printFile is the name of the file where the graph will be printed (.jpg)
%OUTPUTS
    %dist is an array containing 2 columns. 1st column is the gene 
    %intensity and the 2nd column is the reaction flux. 
%Narayanan Sadagopan 10/12/13

dist = zeros(length(intenst),2);
count = 1;

%create temporary file
fileID=fopen(fileName,'r');
c1 = textscan(fileID, '%s %s %s',1);
C = textscan(fileID, '%f %f %f');
hdr={'gene','mean','var'}; 
txt=sprintf('%s\t',hdr{:});
txt(end)='';
dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter','\t','precision',10);
fclose(fileID);

%modify gene intensity and run Falcon
for x = 1:length(intenst)

    %make change and rewrite file
    C{1,2}(1) = intenst(x); 
    dlmwrite('tempFileForAnalyzeFalcon.csv',txt,'');
    dlmwrite('tempFileForAnalyzeFalcon.csv',C,'-append','delimiter','\t','precision',10);
    
    %run Falcon
    [vIrrev vRev] = runFalcon(recMod,'tempFileForAnalyzeFalcon.csv', rc);
    dist(count,1) = intenst(x);
    dist(count,2) = vRev(rxnOfInt);
    count = count + 1;
    
end

delete('tempFileForAnalyzeFalcon.csv');

%plot
figure;
plot(dist(:,1),dist(:,2),'b-o');
title('Flux vs Intensity');
xlabel(sprintf('intensity of the gene %d', C{1,1}(1)));
ylabel(sprintf('flux for reaction %d', rxnOfInt));
grid on;
print('-dpng',printFile);
close(figure);

   
