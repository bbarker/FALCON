function output = makeMethodFigures(figName)

output = [];

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


if strcmp(figName, 'CorrSigma')
    rec2Title = ['Human Recon2 flux sensitivity to gene noise' char(10)];
    %rec2CS = load('pertData_humanrecon2_K562_0.5_1000_2014121142510012932_1000rec203_1.mat')
    rec2CS = load('pertData_humanrecon2_K562_4_1000_20141220363133847_1000_1.mat');
    CorrSigma(rec2CS.sigmaVec, rec2CS.PcorrV, ...
        [rec2Title 'default constraints'], 'Pearson''s r');

    rec2CS_1C = load('pertData_humanrecon2_K562_4_1000_201412314444710728_med_1000_1.mat');
    CorrSigma(rec2CS_1C.sigmaVec, rec2CS_1C.PcorrV, ...
        [rec2Title 'RPMI-Medium constrained'], 'Pearson''s r');

    rec2CS_2C = load('pertData_humanrecon2_K562_4_1000_20141232331358614_med_coreSign_1000_1.mat');
    CorrSigma(rec2CS_2C.sigmaVec, rec2CS_2C.PcorrV, ...
        [rec2Title 'RPMI-Medium and CORE-Sign constrained'], 'Pearson''s r');

    rec2CS_3C = load(['pertData_humanrecon2_K562_4_1000_201412324558226941_' ...
                      'rec203_med_coreSign_imputed_NotAll_1000_1.mat']);
    CorrSigma(rec2CS_3C.sigmaVec, rec2CS_3C.PcorrV, ...
        [rec2Title 'RPMI-Medium and CORE-Imputed Enzymatic-Reaction directionality constrained'], 'Pearson''s r');

    rec2CS_4C = load(['pertData_humanrecon2lmomaall_K562_4_1000_2014124234511045576' ...
                      '_rec203_med_coreSign_imputed_All_1000_1.mat']);
    CorrSigma(rec2CS_4C.sigmaVec, rec2CS_4C.PcorrV, ...
        [rec2Title 'RPMI-Medium and CORE-Imputed All-Reaction directionality constrained'] , 'Pearson''s r');


    y7dCS = load('pertData_yeast700cobra_genedata_75_4_1000_201412112022547807_1000y7dir_1.mat');
    CorrSigma(y7dCS.sigmaVec, y7dCS.PcorrV, ...
        'Yeast 7 (minimally constrained) flux sensitivity to gene noise', 'Pearson''s r');
    y7ndCS = load('pertData_yeast700cobra_genedata_75_4_1000_201412121527771331_100y7noDir_1.mat');
    CorrSigma(y7ndCS.sigmaVec, y7ndCS.PcorrV, ...
        'Yeast 7 (highly constrained) flux sensitivity to gene noise', 'Pearson''s r');

    %Enzyme Complexities:
    CorrSigma(y7dCS.sigmaVec, y7dCS.PcorrE, ...
        'Yeast 7 Enzyme Complex sensitivity to gene noise', 'Pearson''s r');
    CorrSigma(rec2CS.sigmaVec, rec2CS.PcorrE, ...
        'Human Recon2 Enzyme Complex sensitivity to gene noise', 'Pearson''s r');


end

%This needs to be rewritten to work in function form.
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
    legend('Yeast 7', 'Human Recon2');    
end


if strcmp(figName, 'fluxEvalCmpScatter')
    % data comes from compareFluxByEnzymeComplexation.m
    y7ndFCMP = load('FluxByECcomp_yeast_7.00_cobra.xmlgenedata_75.txtnoDir_1000.mat');
    y7dFCMP = load('FluxByECcomp_yeast_7.00_cobra.xmlgenedata_75.txtdir_1000.mat');

    rec2FCMP = load('FluxByECcomp_human_recon_2K562.csvrec203_1000.mat');
    rec2FCMP_3C = load('FluxByECcomp_human_recon_2K562.csvmed_coreSign_imputed_NotAll_1000.mat');

    fluxCmpScatter(y7ndFCMP, 'Yeast7 Highly Constrained', 1/2) ;
    fluxCmpScatter(y7dFCMP, 'Yeast7 Minimally Constrained', 1/2)
    fluxCmpScatter(rec2FCMP, 'Human Recon2 (default constraints)', 1/2, 0, 1);
    fluxCmpScatter(rec2FCMP_3C, 'Human Recon2 (highly constrained)', 1/2, 0, 1);
end

if strcmp(figName, 'fluxGrpCmpScatter')
    % data comes from compareFluxByEnzymeComplexation.m
    
    y7ndFCMP = load('FluxByGroupComp_yeast_7.00_cobra.xmlgenedata_75.txtrgroupTest_y7nd1k_1000.mat');
    y7dFCMP = load('FluxByGroupComp_yeast_7.00_cobra.xmlgenedata_75.txtrgroupTest_y7d1k_1000.mat');

    fluxGrpCmpScatter(y7ndFCMP, 'Yeast7 Highly Constrained', 1/2) ;
    fluxGrpCmpScatter(y7dFCMP, 'Yeast7 Minimally Constrained', 1/2)
end

% data generated from yeastResults.m
% Be sure to set dC = to dataCells75 or dataCells85
if strcmp(figName, 'fluxBarsTables')
    y5d_75  = importdata('genedata_75.txt_results_all_Rep100_y5dir_735619.813.csv');
    y5nd_75 = importdata('genedata_75.txt_results_all_Rep100_y5orig_735604.7928.csv');
    y7d_75  = importdata('genedata_75.txt_results_all_Rep100_y7dir_735604.8936.csv');
    y7nd_75 = importdata('genedata_75.txt_results_all_Rep100_y7orig_735604.8164.csv');

    % y5d_85  = importdata('');
    y5nd_85 = importdata('genedata_85.txt_results_all_Rep100_y5orig_735604.7929.csv');
    y7d_85  = importdata('genedata_85.txt_results_all_Rep100_y7dir_735604.8949.csv');
    y7nd_85 = importdata('genedata_85.txt_results_all_Rep100_y7orig_735604.8169.csv ');

    dataCells75 = {y5d_75, 'Yeast5 MC'; ...
                   y7d_75, 'Yeast7 MC'; ...
                   y5nd_75, 'Yeast5 HC'; ...
                   y7nd_75, 'Yeast7 HC'
    };

    %Need a fancy work around to support multiline xticklabels:
    % http://www.mathworks.com/matlabcentral/answers/101922
    %dtaCells75 = {y5d_75, ['Yeast5'  char(10)  'min. constrained']; ...
    %              y7d_75, ['Yeast7'  char(10) 'min. constrained']; ...
    %              y5nd_75, ['Yeast5' char(10)  'highly constrained']; ...
    %              y7nd_75, ['Yeast7' char(10)  'highly constrained']
    %;

 
    nFlux = 7;
    methCols = [1 3 4 5 6 7 9];
    methNames = {'Experimental', 'Standard FBA', 'Fitted FBA', 'GIMME', 'iMAT', ...
        'Lee et al.', 'FALCON'};
    nMeth = length(methCols);
    %Define data cols: 
    
    hasSTD = [7, 9];
    %create data matrix:

    % What do we want: 
    % Groups = Model + Growth Media (2 * 2)
    % BarColorInGroup = Method + Experimental

    dC = dataCells75;
    ndC = length(dC);
    for i = 1:nFlux
    %for i = 1:1
       metTitle = y5d_75.textdata{i+1, 1};
       dMean = zeros(ndC, nMeth);
       dSTD  = zeros(ndC, nMeth);
       for j = 1:nMeth
           for k = 1:ndC
               methIdx = methCols(j);
               dMean(k, j) = dC{k}.data(i, methIdx);
               if any(methIdx == hasSTD)
                   dSTD(k, j) = dC{k}.data(i, methIdx + 1);
               else
                   dSTD(k, j) = 0;
               end
           end % end for k
       end % end for j
       figure();
       set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
       set(gca, 'FontSize', 23);
       barwitherr(dSTD, [1:ndC], dMean);
       xlabels = dC(:, 2);
       output = xlabels;
       set(gca,'XTickLabel', xlabels);
       ylabel('Flux (mmol/gDW/h');
       title(dC{1}.textdata{i+1});
       legend(methNames);
       colormap gray;
    end % end for i
end

if strcmp(figName, 'modelTime')
    1;
end
end % end of makeMethodFigures



%%%   Convenience functions   %%%

function outliers = getOutliers(vec, nOut, message, otherVecs)
% otherVecs is a list of data to print associated to the outlier
sortedVec = sort(vec(:)');
Obound = sortedVec(end-nOut);
ok = find(vec < Obound);
outliers = find(vec >= Obound);
primaryData = [num2str(outliers(:)') char(10)];
allData = '';
for i = 1:length(otherVecs)
    V = otherVecs{i};
    allData = [allData num2str(V(outliers)') char(10)];
end 
allData = [primaryData allData];
% print outlier information:
disp([message char(10) allData]);
end % end of [getOutliers]

function fluxCmpScatter(ss, T, sRat, nXout, nYout, nXEout, nYEout)
% nOut is the number of outliers to remove
% sRat is fraction of stdev on one side of error bar.
if ~exist('nXout', 'var') nXout = 0; end
if ~exist('nYout', 'var') nYout = 0; end
if ~exist('nXEout', 'var') nXEout = 0; end
if ~exist('nYEout', 'var') nYEout = 0; end
stdOutliers = [];
mXoutliers  = []; mYoutliers  = [];
mXEoutliers = []; mYEoutliers = [];
if nXout > 0
    msg = 'Reaction indices with outlier x-values: ';
    mXoutliers = getOutliers(ss.v_lee, nXout, msg, {ss.v_lee_s});
end
if nYout > 0
    msg = 'Reaction indices with outlier y-values: ';
    mYoutliers = getOutliers(ss.v_md, nYout, msg, {ss.v_md_s});
end
if nXEout > 0
    msg = 'Reaction indices with STD outlier x-values: ';
    mXEoutliers = getOutliers(ss.v_lee_s, nXEout, msg, {ss.v_lee});
end
if nYEout > 0
    msg = 'Reaction indices with STD outlier y-values: ';
    mYEoutliers = getOutliers(ss.v_md_s, nYEout, msg, {ss.v_md});
end

scatterError(ss.v_lee, ss.v_md, sRat*ss.v_lee_s, sRat*ss.v_md_s, ...
    union(mXoutliers, mYoutliers), mXEoutliers, mYEoutliers);
title(T);
xlabel('Flux from direct evaluation');
ylabel('Flux from minimum-disjunction');
end % end of fluxCmpScatter


%Forgive me for code duplication...

function fluxGrpCmpScatter(ss, T, sRat, nXout, nYout, nXEout, nYEout)
% nOut is the number of outliers to remove
% sRat is fraction of stdev on one side of error bar.
if ~exist('nXout', 'var') nXout = 0; end
if ~exist('nYout', 'var') nYout = 0; end
if ~exist('nXEout', 'var') nXEout = 0; end
if ~exist('nYEout', 'var') nYEout = 0; end
stdOutliers = [];
mXoutliers  = []; mYoutliers  = [];
mXEoutliers = []; mYEoutliers = [];
if nXout > 0
    msg = 'Reaction indices with outlier x-values: ';
    mXoutliers = getOutliers(ss.v_md_nogrp, nXout, msg, {ss.v_md_nogrp_s});
end
if nYout > 0
    msg = 'Reaction indices with outlier y-values: ';
    mYoutliers = getOutliers(ss.v_md, nYout, msg, {ss.v_md_s});
end
if nXEout > 0
    msg = 'Reaction indices with STD outlier x-values: ';
    mXEoutliers = getOutliers(ss.v_md_nogrp_s, nXEout, msg, {ss.v_md_nogrp});
end
if nYEout > 0
    msg = 'Reaction indices with STD outlier y-values: ';
    mYEoutliers = getOutliers(ss.v_md_s, nYEout, msg, {ss.v_md});
end

scatterError(ss.v_md_nogrp, ss.v_md, sRat*ss.v_md_nogrp_s, sRat*ss.v_md_s, ...
    union(mXoutliers, mYoutliers), mXEoutliers, mYEoutliers);
title(T);
xlabel('Flux from FALCON: no reaction groups');
ylabel('Flux from FALCON');
end % end of fluxGrpCmpScatter


function CorrSigma(cX, cY, T, ctype)
    [xI, yI] = intervalAggregate(cX, cY, @median, (max(cX)-min(cX))/10, 0.3); 

    figure();
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    set(gca, 'FontSize', 23);
    hold all;     
    scatter(cX, cY);
    plot(xI, yI, 'color', 'g', 'LineWidth', 3);

    title(T);
    xlabel('\sigma (error distribution parameter)');
    ylabel(ctype);
end % end of CorrSigma
