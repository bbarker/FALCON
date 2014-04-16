function output = makeMethodFigures(figName)

output = [];

% Whether to draw titles in figures.
drawTitles = false;

%Need to annotate (at least partly) which figure numbers
%correspond to which figName(s) below.


% permuted experimental correlation distribution
% see randomExpressionSim.m to generate data
if strcmp(figName, 'pExp') % (figure YpermCorr)
    %Be sure to declare these as globals in the main workspace first.
    s75 = load('run1randCorr_5000_genedata_75.mat', '-mat', 'PearsonExpVec');
    s85 = load('run1randCorr_5000_genedata_85.mat', '-mat', 'PearsonExpVec');
    expTitle = '';
    if drawTitles
        expTitle = ['Predictions for actual versus permuted expression:' ...
                     char(10) 'experimental flux correlation'];
    end
    plotSmoothHistPairs(s75.PearsonExpVec, s85.PearsonExpVec, expTitle,  ...
        'Pearson''s r', -1, 1.3, -1,                                     ...
         0.99, '75% Max Growth', 0.98, '85% Max Growth');


    s75 = load('dirr1y7randCorr_5000_genedata_75.mat', '-mat', 'PearsonExpVec');
    s85 = load('dirr1y7randCorr_5000_genedata_85.mat', '-mat', 'PearsonExpVec');
    plotSmoothHistPairs(s75.PearsonExpVec, s85.PearsonExpVec, expTitle,  ...
        'Pearson''s r', -1, 1.3, 0.008,                                  ...
         0.98, '75% Max Growth', 0.97, '85% Max Growth');

end
% permuted flux correlation distribution
% see randomExpressionSim.m to generate data
if strcmp(figName, 'pVkV')

    kw = 0.01;
    expTitle = '';
    if drawTitles
        expTitle = ['Flux vector correlation with fluxes' char(10) ...
                     'derived from permuted expression'];
    end
    s75 = load('run1randCorr_5000_genedata_75.mat', '-mat', ...
               'PearsonVec', 'KendallVec');
    s85 = load('run1randCorr_5000_genedata_85.mat', '-mat', ...
               'PearsonVec', 'KendallVec');
    plotSmoothHistPairs(s75.PearsonVec, s85.PearsonVec, expTitle, ...
        'Pearson''s r', -0.55, 1, kw);
    plotSmoothHistPairs(s75.KendallVec, s85.KendallVec, expTitle, ...
        'Kendall''s tau', -0.55, 1, kw);
    %
    s75 = load('dirr1y7randCorr_5000_genedata_75.mat', '-mat', ...
               'PearsonVec', 'KendallVec');
    s85 = load('dirr1y7randCorr_5000_genedata_85.mat', '-mat', ...
               'PearsonVec', 'KendallVec');
    plotSmoothHistPairs(s75.PearsonVec, s85.PearsonVec, expTitle, ...
        'Pearson''s r', -0.55, 1, kw);
    plotSmoothHistPairs(s75.KendallVec, s85.KendallVec, expTitle, ...
        'Kendall''s tau', -0.55, 1, kw);


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
    rec2Title = ['Human Recon 2 flux sensitivity to gene noise' char(10)];
    %rec2CS = load('pertData_humanrecon2_K562_0.5_1000_2014121142510012932_1000rec203_1.mat')
    rec2CS = load('pertData_humanrecon2_K562_4_1000_20141220363133847_1000_1.mat');
    CorrSigma(rec2CS.sigmaVec, rec2CS.PcorrV, ...  % 1
        [rec2Title 'default constraints'], 'Pearson''s r');

    rec2CS_1C = load('pertData_humanrecon2_K562_4_1000_201412314444710728_med_1000_1.mat');
    CorrSigma(rec2CS_1C.sigmaVec, rec2CS_1C.PcorrV, ... % 2
        [rec2Title 'RPMI-medium constrained'], 'Pearson''s r');

    rec2CS_2C = load('pertData_humanrecon2_K562_4_1000_20141232331358614_med_coreSign_1000_1.mat');
    CorrSigma(rec2CS_2C.sigmaVec, rec2CS_2C.PcorrV, ... % 3
        [rec2Title 'RPMI-medium and CORE-sign constrained'], 'Pearson''s r');

    rec2CS_3C = load(['pertData_humanrecon2_K562_4_1000_201412324558226941_' ...
                      'rec203_med_coreSign_imputed_NotAll_1000_1.mat']);
    CorrSigma(rec2CS_3C.sigmaVec, rec2CS_3C.PcorrV, ... % 4
        [rec2Title 'RPMI-medium and CORE-imputed enzyme-reaction directionality constrained'], 'Pearson''s r');

    rec2CS_4C = load(['pertData_humanrecon2lmomaall_K562_4_1000_2014124234511045576' ...
                      '_rec203_med_coreSign_imputed_All_1000_1.mat']);
    CorrSigma(rec2CS_4C.sigmaVec, rec2CS_4C.PcorrV, ... % 5
        [rec2Title 'RPMI-medium and CORE-imputed all-reaction directionality constrained'] , 'Pearson''s r');


    y7dCS = load('pertData_yeast700cobra_genedata_75_4_1000_201412112022547807_1000y7dir_1.mat');
    CorrSigma(y7dCS.sigmaVec, y7dCS.PcorrV, ... % 6
        'Yeast 7 (minimally constrained) flux sensitivity to gene noise', 'Pearson''s r');
    y7ndCS = load('pertData_yeast700cobra_genedata_75_4_1000_201412121527771331_100y7noDir_1.mat');
    CorrSigma(y7ndCS.sigmaVec, y7ndCS.PcorrV, ... % 7
        'Yeast 7 (highly constrained) flux sensitivity to gene noise', 'Pearson''s r');

    %Enzyme Complexities:
    CorrSigma(y7dCS.sigmaVec, y7dCS.PcorrE, ... % 8
        'Yeast 7 enzyme complex sensitivity to gene noise', 'Pearson''s r');
    CorrSigma(rec2CS.sigmaVec, rec2CS.PcorrE, ... % 9
        'Human Recon 2 enzyme complex sensitivity to gene noise', 'Pearson''s r');


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
    if drawTitles
        title([cType ' correlations for perturbed expression data']);
    end
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

%
% Data from this file is generated by 
%
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
    if drawTitles
        title(['Differences in complex abundance ' char(10) ...
               'between minimum-disjunction and direct evaluation']);
    end
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
    fluxCmpScatter(rec2FCMP, 'Recon2 (default constraints)', 1/2, 0, 1);
    fluxCmpScatter(rec2FCMP_3C, 'Recon2 (highly constrained)', 1/2, 0, 1);
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
    y5d_75_f = importdata('genedata_75.txt_results_mindisj2_Y5_dir100_735654.4775.csv');

    y5nd_75 = importdata('genedata_75.txt_results_all_Rep100_y5orig_735604.7928.csv');
    y5nd_75_f = importdata('genedata_75.txt_results_mindisj2_Y5_NoDir100_735654.4724.csv');

    y7d_75  = importdata('genedata_75.txt_results_all_Rep100_y7dir_735604.8936.csv');
    y7d_75_f = importdata('genedata_75.txt_results_mindisj2_Y7_dir100_735653.8628.csv');

    y7nd_75 = importdata('genedata_75.txt_results_all_Rep100_y7orig_735604.8164.csv');
    y7nd_75_f = importdata('genedata_75.txt_results_mindisj2_Y7_NoDir100_735653.8996.csv');

    % y5d_85  = importdata('');
    y5nd_85 = importdata('genedata_85.txt_results_all_Rep100_y5orig_735604.7929.csv');
    y7d_85  = importdata('genedata_85.txt_results_all_Rep100_y7dir_735604.8949.csv');
    y7nd_85 = importdata('genedata_85.txt_results_all_Rep100_y7orig_735604.8169.csv ');

    dataCells75 = {y5d_75 'Yeast5 MC'; ...
                   y7d_75 'Yeast7 MC'; ...
                   y5nd_75 'Yeast5 HC'; ...
                   y7nd_75 'Yeast7 HC'
    };

    dataCells85 = {... %y5d_85 'Yeast5 MC'; ...
                   y7d_85 'Yeast7 MC'; ...
                   y5nd_85 'Yeast5 HC'; ...
                   y7nd_85 'Yeast7 HC'
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
   
    %use updated mindisj data:
    y5d_75.data(:, 9:10)  = y5d_75_f.data(:, 9:10);
    y5nd_75.data(:, 9:10) = y5nd_75_f.data(:, 9:10);
    y7d_75.data(:, 9:10)  = y7d_75_f.data(:, 9:10); 
    y7nd_75.data(:, 9:10) = y7nd_75_f.data(:, 9:10);

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
       set(gca, 'FontSize', 28);
       barwitherr(sgnLog10p1(dSTD), [1:ndC], sgnLog10p1(dMean));
       xlabels = dC(:, 2);
       set(gca,'XTickLabel', xlabels);
       %ylabel('Log10(1 + Flux (mmol/gDW/h))');
       ylabel('Flux');
       if drawTitles
           title(dC{1}.textdata{i+1});
       end
       colormap(gray);
    end % end for i
    % One for a legend:
    figure();
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    set(gca, 'FontSize', 28);
    barwitherr(dSTD, [1:ndC], dMean);
    legend(methNames);
    colormap(gray)

    % Now we make a table for the average timing and correlation
    % cols: number of methods + 2, rows: 1 + length(dC)
    dC = [[dataCells75 repmat({'75% max growth'}, length(dataCells75), 1)]; ...
          [dataCells85 repmat({'85% max growth'}, length(dataCells85), 1)]];
    ndC = length(dC);
    pearsonTab = cell(ndC + 1, nMeth + 2); 
    timingTab  = cell(ndC + 1, nMeth + 2); 
    pearsonTab{1, 1} = 'Max Growth Rate Percent';
    timingTab{1, 1}  = 'Max Growth Rate Percent';
    pearsonTab{1, 2} = 'Model';
    timingTab{1, 2}  = 'Model';
    sz_pTab = size(pearsonTab)
    sz_dC = size(dC)
    pearsonTab(2:end, 1) = dC(:, 3);
    timingTab(2:end, 1)  = dC(:, 3);
    pearsonTab(2:end, 2) = dC(:, 2);
    timingTab(2:end, 2)  = dC(:, 2);
    pearsonTab(1, 3:end) = methNames;
    timingTab(1, 3:end) = methNames;

    for j = 1:nMeth
        methIdx = methCols(j);
        for k = 1:ndC
            pearsonTab{1 + k, 2 + j} = num2str(dC{k}.data(8, methIdx)); 
            timingTab{1 + k, 2 + j}  = num2str(dC{k}.data(9, methIdx)); 
        end
    end % for j
    output.pearsonTab = pearsonTab;
    output.timingTab = timingTab;
end

if strcmp(figName, 'modelTime')
    1;
end



%%%   Convenience functions   %%%

function xt = sgnLog10p1(x)
    lgtrans = @(y) sign(y) * log10(1+abs(y));
    xt = arrayfun(lgtrans, x);
end % of sgnLog10p1

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
if drawTitles
    title(T);
end
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
if drawTitles
    title(T);
end
xlabel('Flux from FALCON: no reaction groups');
ylabel('Flux from FALCON');
end % end of fluxGrpCmpScatter


function CorrSigma(cX, cY, T, ctype)
    [xI, yI] = intervalAggregate(cX, cY, @median, (max(cX)-min(cX))/10, 0.3); 

    figure();
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    set(gca, 'FontSize', 32);
    hold all;     
    scatter(cX, cY);
    plot(xI, yI, 'color', 'g', 'LineWidth', 3);
    if drawTitles
        title(T);
    end
    xlabel('\sigma (error distribution parameter)');
    ylabel(ctype);
end % end of CorrSigma

end % end of makeMethodFigures
