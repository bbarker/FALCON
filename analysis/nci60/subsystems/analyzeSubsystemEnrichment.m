%Takes a particular subsystem in human recon 2 and converts Entrez gene ids
%to IPI ids. Looks at quantified protein and mRNA amounts in data set
%for the IPI ids.One potential flaw, and/or statements in gene rules 
function [num num2 pMRNA pProtein missingGenesMRNA missingGenesProtein residuals]=analyzeSubsystemEnrichment(converter, data, subsys, rec2)

%INPUT
%converter is a (1x2) cell array with Entrez ids in {1,1} and the corresponding IPI ids in {1,2}.
%data is the output of readModelExpression
%subsys is a human recon 2 subsystem
%rec2 is the human recon 2 model
%OUTPUT
%num is number of mRNA intensities
%num2 is number of protein intensities
%pMRNA is genes in paper with mRNA intensity/possible genes for the subsystem
%pProtein is genes in paper with protein intensity/possible genes for the subsystem
%missingGenesMRNA a cell array containing the Entrez id of all genes in the 
%model for the subsystem missing an MRNA intensity
%missingGenesProtein a cell array containing the Entrez id of all genes in 
%the model for the subsystem missing an protein intensity
%residual is vector containing the absolute value of the difference between
%MRNA and protein intensites
%Narayanan Sadagopan, September 2013

num = 0;   %number of elements of mRNA
num2 = 0;  %number of elements of protein
residuals(1) = 0; %difference between protein and mRNA intensity values
allData(1) = 5; %all mRNA intensities    
allData2(1) = 5;  %all protein intensities                
pMRNA = 0;  %genes in paper with mRNA intensity/possible genes for that subsystem
pProtein = 0;  %genes in paper with protein intensity/possible genes for that subsystem
missingGenesMRNA{1} = '0';  %genes in model that don't have an mRNA intensity
missingGenesProtein{1} = '0';  %genes in model that don't have an protein intensity
               
count = 1;  %indexing for mRNA intensities
count2 = 1;  %indexing for protein intensities 
count3 = 1;  %indexing for residuals
count4 = 0;  %counts total possible genes for the subsystem
count5 = 0;  %counts total genes with an mRNA intensity
count6 = 0;  %counts total genes with an protein intensity
count9 = 1;  %indexing for missing mRNA genes
count10 = 1;  %indexing for missing protein genes

for x = 1:7440
    if (strcmp(rec2.subSystems{x}, subsys))
        for b = 1:2194
            if (rec2.rxnGeneMat(x,b) == 1)
                count4 = count4 + 1;    
                truth3 = 0; %to check for gene's presence in mRNA
                truth4 = 0; %to check for gene's presence in protein   
                for y = 1:length(converter{1,1})
                    if (strcmp(converter{1,1}(y),rec2.genes(b)))
                        for a = 1:length(data{1,2})
                            truth = 0; %to check for mRNA intensity for residual calc
                            truth2 = 0; %to check for protein intensity for residual calc   
                            if (strcmp(converter{1,2}{y}, data{1,2}{a}))
                                if (data{1,5}(a) == data{1,5}(a))  %removes any NaN values
                                    allData(count) = data{1,5}(a);
                                    count = count + 1;
                                    truth = 1;
                                    truth3 = 1;
                                end
                                if (data{1,6}(a) == data{1,6}(a))  %removes any NaN values
                                    allData2(count2) = data{1,6}(a);
                                    count2 = count2 + 1;
                                    truth2 = 1;
                                    truth4 = 1;
                                end
                            end
                            if (truth&&truth2)
                                residuals(count3) = abs(allData(count-1) - allData2(count2-1));
                                count3 = count3 + 1;
                            end
                        end
                    end
                end
                if (~truth3)
                    missingGenesMRNA{count9} = rec2.genes(b);
                    count9 = count9 + 1;
                else
                    count5 = count5 + 1;
                end
                if (~truth4)
                    missingGenesProtein{count10} = rec2.genes(b);
                    count10 = count10 + 1;
                else
                    count6 = count6 + 1; 
                end
            end
        end
    end
end
num = length(allData);
if(allData(1) == 5) %if no mRNA intensites are present, length set to 0   
    num = 0;
end
num2 = length(allData2);
if(allData2(1) == 5) %if no protein intensites are present, length set to 0  
    num2 = 0;
end
pMRNA = count5./count4;
pProtein = count6./count4;
if (count4 == 0) %to prevent nan values
    pMRNA = 0;
    pProtein = 0;
end
