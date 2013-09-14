%generate a histogram for protein and mRNA percent gene coverage over all subsystems 
function [] = subsystemsMakeHist(converter,data,subsys,rec2,name,name2)

%INPUT
%converter is a (1x2) cell array with Entrez ids in {1,1} and the corresponding IPI ids in {1,2}.
%data is the output of readModelExpression
%subsys is ALL human recon 2 subsystems
%rec2 is the human recon 2 model
%name is the filename where mRNA histogram is generated (ex. 'hist.jpg');
%name2 is the filename where protein histogram is generated (ex. 'hist.jpg');
%OUTPUT
%generates two files containing histograms
%Narayanan Sadagopan, September 2013

%get data of gene coverage for all subsystems
count = 1;
for x = 1:length(subsys)
    [num{count} num2{count} pMRNA{count} pProtein{count}]=analyzeSubsystemEnrichment(converter,data,subsys{x},rec2);
    count = count + 1;
end

%switch format for histogram
m = [pMRNA{:}];
p = [pProtein{:}];

%make histograms
fh = figure;
hist(m,20);
title_handle = title('MRNA');
set(title_handle,'String','MRNA');
print(fh,'-dpng',name);
fclose(fh);
fh2 = figure;
hist(p,20);
title_handle = title('Protein');
set(title_handle,'String','Protein');
print(fh2,'-dpng',name2);
fclose(fh2);

