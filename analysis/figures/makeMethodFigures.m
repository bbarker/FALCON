%Need to annotate (at least partly) which figure numbers
%correspond to which figName(s) below.


% permuted experimental correlation distribution
% see randomExpressionSim.m to generate data
if strcmp(figName, 'pExp')
    plotSmoothHistPairs(PearsonExpVec75, PearsonExpVec85,              ...
        ['Predictions for actual versus permuted expression:' char(10) ...
        'experimental flux correlation'] ,  ...                         ...
                'Pearson''s r', 0.9820, '75% Max Growth', 0.9527, '85% Max Growth')

%        'Pearson''s r', 0.9879, '75% Max Growth', 0.9664, '85% Max Growth') 
% unconstrained: 
end
% permuted flux correlation distribution
% see randomExpressionSim.m to generate data
if strcmp(figName, 'kV')
    plotSmoothHistPairs(KendallVec75, KendallVec85,        ...
        ['Flux vector correlation with fluxes' char(10)    ...
        'derived from permuted expression'], 'Kendall''s tau')
end
if strcmp(figName,'pV')
    plotSmoothHistPairs(PearsonVec75, PearsonVec85,        ...
        ['Flux vector correlation with fluxes' char(10)    ...
        'derived from permuted expression'], 'Pearson''s r')
end

% create expression/flux comparison figures.
% see compareExpressionToFlux.m to generate data matrix
if strcmp(figName, 'EvVHum')
    %load colormap
    load('expFluxCmap.mat');
    fig = figure;
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    imagesc(corrMatK562overMDAMB231);
    colormap(expFluxCmap);    
    colorbar;    
    set(gca, 'FontSize', 23);
    set(gca,'Xtick',1:7,'XTickLabel',{'|Flux|', 'EC', ...
        'sum', 'mean', 'median', 'min', 'max'});    
    set(gca,'Ytick',1:7,'YTickLabel',{'|Flux|', 'EC', ...
        'sum', 'mean', 'median', 'min', 'max'});
    axis square;
end
if strcmp(figName, 'EvVY')
    %load colormap; identical to above except diff mat - make a fun?
    load('expFluxCmap.mat');
    fig = figure;
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    imagesc(corrMat75over85);
    colormap(expFluxCmap);    
    colorbar;    
    set(gca, 'FontSize', 23);
    set(gca,'Xtick',1:7,'XTickLabel',{'|Flux|', 'EC', ...
        'sum', 'mean', 'median', 'min', 'max'});    
    set(gca,'Ytick',1:7,'YTickLabel',{'|Flux|', 'EC', ...
        'sum', 'mean', 'median', 'min', 'max'});
    axis square;    
end

if strcmp(figName, 'PertAn')
    %e.g.:
    %corrE = PcorrE; corrV=PcorrV; cType = 'Pearson'; figName = 'PertAn'; makeMethodFigures
    fig = figure;
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    hold all;    
    scatter(corrE, corrV);
    set(gca, 'FontSize', 32);
    xlim([-0.1 1.001]);
    hls = lsline;
    set(hls, 'LineWidth', 3);
    title([cType ' correlations for perturbed expression data']);
    xlabel([cType '''s correlation between expression'] );
    ylabel([cType '''s correlation between fluxes'] );
    mb = polyfit(corrE, corrV, 1);
    r = corr(corrE', corrV');
    xA = [0.55 0.55];
    %yA = [-0.2 -0.2]; %YN
    yA = [-0.1 -0.1]; %rec2
    [xF, yF] = ds2nfu(xA, yA);    
    a1 = annotation('textbox', [xF(1) yF(1) .1 .1], 'FontSize', 24, 'String', ...
        ['y = ' num2str(mb(1)) 'x + ' num2str(mb(2)) ', r = ' num2str(r)]);        
end

if strcmp(figName, 'expCompare')
    iMax = 13;
    yDat = load('compareEnzymeExpression_y7_default1000.mat');
    hDat = load('compareEnzymeExpression_hrec2_default200.mat');
    fig = figure;
    set(gcf, 'Position', get(0,'Screensize'));
    plot(1:iMax, 100*yDat.nnanDiffTotal(1:iMax)/ yDat.nnanTotal, 'LineWidth', 4);
    set(gca, 'FontSize', 26);
    xlim([1 iMax]);
    set(gca,'XTick',[1:2:iMax]);
    title(['Differences in complex abundance ' char(10) ...
           'between minimum-disjunction and direct evaluation']);
    xlabel('Number of permuted expression vectors');
    ylabel('Percentage of ezymatic reactions');
    hold all;
    plot(1:iMax, 100*hDat.nnanDiffTotal(1:iMax)/ hDat.nnanTotal, 'LineWidth', 4);
    legend('Yeast 7', 'Human Recon 2');    
end

if strcmp(figName, 'fluxCmpScatter')


end

if strcmp(figName, 'modelTime')
    1;
end
