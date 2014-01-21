function plotSmoothHistPairs(Data1, Data2, figTitle, xLabel, ...
    Est1, Label1, Est2, Label2)

% Requires ds2nfu: 
% http://www.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion

color1 = 'b';
color2 = 'g';
lWidth = 3;

fig = figure;
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
hold all;
set(gca, 'FontSize', 36);

%default:
%[y1, x1] = ksdensity(Data1);
[y1, x1] = ksdensity(Data1, 'width', 0.006 );
plot(x1, y1, 'color', color1, 'LineWidth', lWidth);
%default:
%[y2, x2] = ksdensity(Data2);
[y2, x2] = ksdensity(Data2, 'width', 0.006 );
plot(x2, y2, 'color', color2, 'LineWidth', lWidth);

%default:
%xlim([-1 1.25]);

% Kendall/flux
xlim([-0.19 1.001]);

xlabel(xLabel);
ylabel('Density');
title(figTitle);

if exist('Est1', 'var')
    xA = [Est1 Est1];
    yA = [0.3 0.0];
    [xF, yF] = ds2nfu(xA, yA);
    a1 = annotation('textarrow', xF, yF, ...
        'String', [Label1 char(10) 'Predicted'], 'color', color1, ...
        'LineWidth', lWidth - 1, 'LineStyle', '--', ...
        'fontSize', 32, 'fontWeight', 'bold');
end

if exist('Est2', 'var')
    xA = [Est2 Est2];
    yA = [0.75 0.0];
    [xF, yF] = ds2nfu(xA, yA);
    a1 = annotation('textarrow', xF, yF, ...
        'String', [Label2 char(10) 'Predicted'], 'color', color2, ...
        'LineWidth', lWidth - 1, 'LineStyle', '--', ...
        'fontSize', 32, 'fontWeight', 'bold');
end
