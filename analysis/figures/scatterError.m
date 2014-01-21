function scatterError(x, y, xe, ye)
%Add optional function to transform xe, ye?

% blue circles with red & grey error bars?

nD = length(x);
dotColor = [1 0.3 0.3]; % conservative pink
yeColor = [0, 0.4, 0.8]; % bright navy blue
xeColor = [0.35, 0.35, 0.35]; % not-too-dark grey

dotSize = 10;

figure();
hold all;


for i = 1:nD
    plot([(x(i) - xe(i)) (x(i) + xe(i))], [y(i) y(i)], 'Color', xeColor);
    plot([x(i) x(i)], [(y(i) - ye(i)) (y(i) + ye(i))], 'Color', yeColor);
end

scatter(x, y, dotSize*ones(1, nD), repmat(dotColor, nD, 1));
