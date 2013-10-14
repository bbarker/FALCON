%Takes a particular subsystem in human recon 2. Looks at quantified 
%protein and mRNA amounts in data set based on Entrez Ids. 
%One potential flaw, and/or statements in gene rules 
function [pMRNA pProtein missingGenesMRNA missingGenesProtein totGenes]=analyzeSubsystemEnrichment(data, subsys, rec2, includeZero)

%INPUT
%data is the output of readModelExpression
%subsys is a human recon 2 subsystem
%rec2 is the human recon 2 model
%includeZero is whether or not zero intensities are counted. 1 is for yes. 0 is for no.
%OUTPUT
%pMRNA is genes in paper with mRNA intensity/possible genes for the subsystem
%pProtein is genes in paper with protein intensity/possible genes for the subsystem
%missingGenesMRNA a cell array containing the Entrez id of all genes in the 
%model for the subsystem missing an MRNA intensity
%missingGenesProtein a cell array containing the Entrez id of all genes in 
%the model for the subsystem missing an protein intensity
%totGenes is the total number of genes for the subsystem
%Narayanan Sadagopan, September 2013

pMRNA = 0;  %genes in paper with mRNA intensity/possible genes for that subsystem
pProtein = 0;  %genes in paper with protein intensity/possible genes for that subsystem
missingGenesMRNA{1} = '0';  %genes in model that don't have an mRNA intensity
missingGenesProtein{1} = '0';  %genes in model that don't have an protein intensity
allGenes(1) = 50; %all genes in model for the subsystem
foundGenesMRNA(1) = 3000;  %genes that have mRNA intensity
foundGenesProtein(1) = 3000;  %genes that have protein intensity
totGenes=0; %total number of genes for subsystem

count = 1;  %indexing for mRNA genes
count2 = 1;  %indexing for protein genes
count3 = 1;  %indexing for residuals
count4 = 1;  %indexing for possible genes 
count5 = 1;  %indexing for missing mRNA genes
count6 = 1;  %indexing for missing protein genes

for x = 1:7440
    if (strcmp(rec2.subSystems{x}, subsys))
        for b = 1:2194
            if (rec2.rxnGeneMat(x,b) == 1)
                allGenes(count4) = b;					
                count4 = count4 + 1;    						
                truth3 = 0; %to check for gene's presence in mRNA
                truth4 = 0; %to check for gene's presence in protein   
                for a = 1:length(data{1,1})
                    if (includeZero) 
                        if (strcmp(rec2.genes{b}, data{1,1}{a}))
                            if (data{1,5}(a) == data{1,5}(a))  %removes any NaN values
                                truth3 = 1;
                                foundGenesMRNA(count) = b;
                                count = count + 1;
                            end
                            if (data{1,6}(a) == data{1,6}(a))  %removes any NaN values
                                truth4 = 1;
                                foundGenesProtein(count2) = b;
                                count2 = count2 + 1;
                            end
                        end
                    else 
                        if (strcmp(rec2.genes{b}, data{1,1}{a}))
                            if (data{1,5}(a) == data{1,5}(a) && data{1,5}(a)~=0)  %removes any NaN values
                                truth3 = 1;
                                foundGenesMRNA(count) = b;
                                count = count + 1;
                            end
                            if (data{1,6}(a) == data{1,6}(a) && data{1,6}(a)~=0)  %removes any NaN values
                                truth4 = 1;
                                foundGenesProtein(count2) = b;
                                count2 = count2 + 1;
                            end
                        end
                    end
                end
                if (~truth3)
                    missingGenesMRNA{count5} = rec2.genes(b);
                    count5 = count5 + 1;
                end
                if (~truth4)
                    missingGenesProtein{count6} = rec2.genes(b);
                    count6 = count6 + 1;
                end
            end
        end
    end
end

pMRNA = length(unique(foundGenesMRNA)) / length(unique(allGenes));
pProtein = length(unique(foundGenesProtein)) / length(unique(allGenes));
totGenes=length(unique(allGenes));
if (allGenes(1) == 50) %to prevent nan values
    pMRNA = 0;
    pProtein = 0;
    totGenes=0;
end
if (foundGenesMRNA(1) == 3000) %correction for foundGenesMRNA declaration
    pMRNA = 0;
end
if (foundGenesProtein(1) == 3000) %correction for foundGenesProtein declaration
    pProtein = 0;
end