function makeBarGraph(x,y,xTicks,xTickLabels,yLims,aTitle,saveLoc)
figure('Position',[300, 300, 2000, 500],'Visible','off');
a=bar(x,y,0.5);
set(gca,'XTickLabel',xTickLabels);
set(gca,'XTick',xTicks);
set(gca,'Ylim',yLims);
set(gca,'FontSize',8);
title(aTitle);
xticklabel_rotate;
saveas(a,saveLoc);
end