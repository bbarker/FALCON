function scatterError(x, y, xe, ye, outM, outE, varargin)
%
% Brandon Barker 01/20/2014
%
% outM and outE are predetermined outliers for (x, y) and 
% (xe, ye), respectively.
%
nD = length(x);

%Make these defaults later:
dotColor = [1 0.3 0.3]; % conservative pink
yeColor = [0, 0.4, 0.8]; % bright navy blue
xeColor = [0.35, 0.35, 0.35]; % not-too-dark grey
dotOutEColor = [0 0 0]; % For plotting points without error bars
 
dotSize = 23;

figure();
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gca, 'FontSize', 23);
hold all;

for i = 1:nD
    if numel(find(outE==i)) == 0
        plot([(x(i) - xe(i)) (x(i) + xe(i))], [y(i) y(i)], 'Color', xeColor);
        plot([x(i) x(i)], [(y(i) - ye(i)) (y(i) + ye(i))], 'Color', yeColor);
    end
end

xOK = x;
xOK(union(outM, outE)) = [];
yOK = y;
yOK(union(outM, outE)) = [];
scatter(xOK, yOK, dotSize, dotColor);
%finish this, make solid
scatter(x(outE), y(outE), 1.5*dotSize, dotOutEColor, 'fill');
set(gca, varargin{:});
axis square;
