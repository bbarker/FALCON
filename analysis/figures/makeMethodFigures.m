% permuted experimental correlation distribution
if strcmp(figName, 'pExp')
	plotSmoothHistPairs(PearsonExpVec75, PearsonExpVec85,              ...
		['Predictions for actual versus permuted expression:' char(10) ...
		'experimental flux correlation'] ,                             ...
		'Pearson''s r', 0.9879, '75% Max Growth', 0.9664, '85% Max Growth')
end

% permuted flux correlation distribution
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