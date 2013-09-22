function subsystemsMakeHist(converter, data, subsys, rec2, includeZero, name, name2)
% generate a histogram for protein and mRNA percent gene coverage over all subsystems 
%INPUT
% converter is a (1x2) cell array with Entrez ids in {1,1} and the corresponding IPI ids in {1,2}.
%     (output of getConverter)
% data is the output of readModelExpression
% subsys is ALL human recon 2 subsystems
% rec2 is the human recon 2 model
% includeZero is whether or not zero intensities are counted. 1 is for yes. 0 is for no.
% name is the filename where mRNA histogram is generated (ex. 'hist.jpg');
% name2 is the filename where protein histogram is generated 
%OUTPUT
% generates two files containing histograms
%
% Narayanan Sadagopan, September 2013



%get data of gene coverage for all subsystems
count = 1;
for x = 1:length(subsys)
    [pMRNA(count) pProtein(count)]=analyzeSubsystemEnrichment(converter, data, subsys{x}, rec2, includeZero);
    count = count + 1;
end



%make histograms
figure;
hist(pMRNA,20);
title('mRNA');
xlabel('ratio of genes present in paper to subsystem genes');
ylabel('frequency');
print('-dpng', name);
close;
figure;
hist(pProtein,20);
title('Protein');
xlabel('ratio of genes present in paper to subsystem genes');
ylabel('frequency');
print('-dpng', name2);
close;

