%For a particular subystem, looks at all the reactions.
%For each reaction, looks at the products. Then the function
%finds total # of reactions in model who's reactants are that product.
%Splits the # of reactions per metabolite into 4 regions and plots
%a pie chart. 
function dist = metOverlap (rec2, subsys, pieName)
%INPUTS
% rec2 is the human recon 2 model
% subsys is a human recon 2 subsystem
% pieName is the name of the file to save the pie chart
%OUTPUTS
% dist is a vector containing the number of reactions that have each 
%metabolite

dist(1)= 0;
count = 1;
countSmall = 0;
countMedium = 0;
countLarge = 0;
countVeryLarge = 0;

%for all reactions in the subsystem, find out many reactions
%their metabolites participate in
for y = 1:length(rec2.rxns)
    if (strcmp(rec2.subSystems{y},subsys))
        %looks only at product metabolites
        temp = find(rec2.S(:,y)==1); 
        for x = 1:length(temp)
            %looks for when those metabolites are reactants of other rxns
            dist(count) = length (find(rec2.S(temp(x),:)==-1)); 
            count = count + 1;
        end
    end
end

%create regions for chart
for x = 1:length(dist)
    if (dist(x)<=5)
        countSmall = countSmall + 1;
    elseif (dist(x)>5 && dist(x)<=15)
        countMedium = countMedium + 1;
    elseif (dist(x)>15 && dist(x)<=50)
        countLarge = countLarge + 1;
    else
        countVeryLarge = countVeryLarge + 1;
    end
end

%make pie chart
pieDist = [countSmall countMedium countLarge countVeryLarge];
labels = {'<=5 rxns ';'6-15  ';'16-50  ';'>50  '};
figure;
h = pie(pieDist);
hText = findobj(h,'Type','text');
percentValues = get(hText,'String');
customStrings = strcat(labels,percentValues);
set (hText,{'String'},customStrings);
title(sprintf('Percentage of the metabolites in %s in different # of reactions',subsys));
print('-dpng', pieName);
close;
            