%Narayanan Sadagopan, September 2013
%Takes a particular subsystem in human recon 2 Converts Entrez to IPI. Looks at quantified protein 
%amounts in data set for the corresponding IPI ids.
%Converter is a cell array with Entrez ids in {1,1} and the corresponding IPI ids in {1,2}.
%data comes from running readModelExpression
%subsys is a hr2 subsystem
%rec2 in hr2 model

%function [num avg sd num2 avg2 sd2 residuals]=findSubsystemsInData2(converter,data,subsys,rec2)
%function [num num2 percentMRNA percentProtein residuals avgResid]=findSubsystemsInData2(converter,data,subsys,rec2)
function [num num2 percentMRNA percentProtein missingGenesMRNA missingGenesProtein]=analyzeSubsystemEnrichment(converter,data,subsys,rec2)
avg=0;  %avg intensity of all proteins
sd=0;   %standard deviation
num=0;   %number of elements
avg2=0;
sd2=0;
num2=0;
residuals(1)=0;
allData2(1)=5;				    %generate histogram to file using percentMRNA and percentProtein
allData(1)=5; %all protein intensities       
percentMRNA=0;
percentProtein=0;
avgResid=5;               
count=1;
count2=1;
count3=1;
count4=0;
count5=0;
count6=0;
missingGenesMRNA{1}=0;
missingGenesProtein{1}=0;
count9=1;
count10=1;
for x=1:7440
	if (strcmp(rec2.subSystems{x},subsys))
		for b=1:2194
			if (rec2.rxnGeneMat(x,b)==1)
				count4=count4+1;	
				truth3=0;
				truth4=0;	
				for y=1:length(converter{1,1})
					if (strcmp(converter{1,1}(y),rec2.genes(b)))
						count7=0;
						count8=0;
                    				for a=1:length(data{1,2})
							%truth=0;
							%truth2=0;  
                        				if (strcmp(converter{1,2}{y},data{1,2}{a}))
		        					if (data{1,5}(a)==data{1,5}(a))		      %removes any NaN values
                                					allData(count)=data{1,5}(a);
                                					count=count+1;
									%truth=1;
									count7=1;
									truth3=1;
			 					end
								if (data{1,6}(a)==data{1,6}(a))		      %removes any NaN values
                                					allData2(count2)=data{1,6}(a);
                                					count2=count2+1;
									%truth2=1;
									count8=1;
									truth4=1;
			 					end
							end
							%if (truth&&truth2)
							%	residuals(count3)=abs(allData(count-1)-allData2(count2-1));
							%	count3=count3+1;
							%end
						end
						if (count7)
							count5=count5+1;
						end
						if (count8)
							count6=count6+1;
						end
       		              end
             			end
				if (~truth3)
					missingGenesMRNA{count9}=rec2.genes(b);
					count9=count9+1;
				end
				if (~truth4)
					missingGenesProtein{count10}=rec2.genes(b);
					count10=count10+1;
				end
			end
		end
	end
end
num=length(allData);
%avg=mean(allData);
%sd=std(allData);
if(allData(1)==5)
	num=0;
end
num2=length(allData2);
%avg2=mean(allData2);
%sd2=std(allData2);    
if(allData2(1)==5)
	num2=0;
end
percentMRNA=count5./count4;
percentProtein=count6./count4;
%avgResid=mean(residuals);

%One potential flaw, and/or statements in gene rules



