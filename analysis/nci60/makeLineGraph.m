function makeLineGraph(x,yArray,xTicks,xTickLabels,yLims,aTitle,saveLoc)
figure('Position',[300, 300, 2000, 500],'Visible','off');
allColors={'y','m','c','r','g','b','w','k'};
colors=allColors(4:6);
for i=1:length(yArray)
    plot(x,yArray{i},colors{i});
    hold on;
end
set(gca,'XTickLabel',xTickLabels);
set(gca,'XTick',xTicks);
set(gca,'Ylim',yLims);
set(gca,'FontSize',8);
title(aTitle);
xticklabel_rotate;
saveas(gcf,saveLoc);
end