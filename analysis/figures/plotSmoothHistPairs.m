function plotSmoothHistPairs(Data1, Data2, figTitle, xLabel, ...
    xmin, xmax, ksweight, Est1, Label1, Est2, Label2)

% Requires ds2nfu: 
% http://www.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion

useDefaultKSweight = false;
if  ksweight < 0
    useDefaultKSweight = true;
end
if ~exist('ksweight', 'var')
    useDefaultKSweight = true;
end

color1 = 'b';
color2 = 'g';
lWidth = 5;

fig = figure;
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold all;
set(gca, 'FontSize', 36);

%default:
%[y1, x1] = ksdensity(Data1);
if useDefaultKSweight
    [y1, x1] = ksdensity(Data1);
else
    [y1, x1] = ksdensity(Data1, 'width', ksweight);
end
plot(x1, y1, 'color', color1, 'LineWidth', lWidth);
%default:
%[y2, x2] = ksdensity(Data2);
if useDefaultKSweight
    [y2, x2] = ksdensity(Data2);
else
    [y2, x2] = ksdensity(Data2, 'width', ksweight, 'kernel', 'triangle');
end
plot(x2, y2, 'color', color2, 'LineWidth', lWidth);

%default:
%xlim([-1 1.25]);

% Kendall/flux
%xlim([-0.19 1.001]);

xlim([xmin, xmax]);

xlabel(xLabel);
ylabel('Density');
if length(figTitle) > 0
   title(figTitle);
end

if exist('Est1', 'var')
    xA = [Est1 Est1];
    testval = max([y1 y2])
    yA = [0.4 0.0]/2.5 * max([y1 y2]);
    [xF, yF] = ds2nfu(xA, yA);
    a1 = annotation('textarrow', xF, yF, ...
        'String', [Label1 char(10) 'Predicted'], 'color', color1, ...
        'LineWidth', lWidth - 1, 'LineStyle', '--', ...
        'fontSize', 30, 'fontWeight', 'bold');
end

if exist('Est2', 'var')
    xA = [Est2 Est2];
    yA = [1 0.0]/2.5* max([y1 y2]);
    [xF, yF] = ds2nfu(xA, yA);
    a1 = annotation('textarrow', xF, yF, ...
        'String', [Label2 char(10) 'Predicted'], 'color', color2, ...
        'LineWidth', lWidth - 1, 'LineStyle', '--', ...
        'fontSize', 30, 'fontWeight', 'bold');
end
