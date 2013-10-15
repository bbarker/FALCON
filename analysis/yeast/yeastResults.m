function yeastResults(model, methodList)

% kieran: 26 apr 12
% Brandon Barker 01/19/13 

%model originally used from paper:
%model_name      = 'yeast_5.21_MCISB.xml';
gene_to_scale   = 'glucose transport';
flux_to_scale   = 'D-glucose exchange';
%model           = readCbModel(model_name);

% This seems to have some problems:
% rsquared = @(f,y)(1 - sum((y-f).^2)/sum((y-mean(y)).^2));

mycorr = @(f,y)(corr(f, y, 'type', 'Pearson'));

%% 75%

[reaction_name, experimental, p_gene_exp, p_standard_fba,                  ...
    p_standard_fba_best, p_gimme, p_shlomi, p_fix, p_falcon, timing] =     ...
    yeastAnalysis(model, 'genedata_75.txt', 'experimental_fluxes_75.txt',  ...
        gene_to_scale, flux_to_scale, methodList);

% display
fresults = fopen('results.csv','w');
fprintf(fresults, ['%s' repmat('\t%s', 1, 8) '\n'],               ...
    'Reaction', 'Experimental 75%', 'Expression', 'Standard FBA', ...
    'Fitted FBA', 'Gimme', 'iMAT', 'Expression (FVA Fixed)', 'falcon');
for k = 1:size(reaction_name,1)
    fprintf(fresults, ['%s' repmat('\t%g', 1, 8) '\n'],                        ...
        reaction_name{k,1}, experimental(k), p_gene_exp(k), p_standard_fba(k), ...
        p_standard_fba_best(k), p_gimme(k), p_shlomi(k), p_fix(k), p_falcon(k));
end

fprintf(fresults, ['%s' repmat('\t%g', 1, 8) '\n'], 'R2', 1,                      ...
    mycorr(p_gene_exp, experimental), mycorr(p_standard_fba, experimental),   ...
    mycorr(p_standard_fba_best,experimental), mycorr(p_gimme, experimental),  ...
    mycorr(p_shlomi, experimental), mycorr(p_fix, experimental),              ...
    mycorr(p_falcon, experimental));

fprintf(fresults, ['%s' repmat('\t%g', 1, 8) '\n'], 'Time', 0,       ...
    timing.gene_exp, timing.standard_fba, timing.standard_fba_best,  ...
    timing.gimme, timing.shlomi, timing.fix, timing.falcon);


return;

%% 85%

disp(' ')

[reaction_name,experimental,p_gene_exp,p_standard_fba,p_standard_fba_best,p_gimme,p_shlomi,p_fix,p_falcon] = ...
    analysis(model, 'genedata_85.txt','experimental_fluxes_85.txt',gene_to_scale,flux_to_scale);

% display
fprintf(fresults,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Reaction','Experimental 85%','Expression','Standard FBA','Fitted FBA','Gimme','iMAT', 'Expression (FVA Fixed)','falcon');
for k = 1:size(reaction_name,1)
    fprintf(fresults,'%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',reaction_name{k,1},experimental(k),p_gene_exp(k),p_standard_fba(k),p_standard_fba_best(k),p_gimme(k),p_shlomi(k),p_fix(k),p_falcon(k));
end
fprintf(fresults,'%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n','R2',1,mycorr(p_gene_exp,experimental),mycorr(p_standard_fba,experimental),mycorr(p_standard_fba_best,experimental),mycorr(p_gimme,experimental),mycorr(p_shlomi,experimental),mycorr(p_fix,experimental),mycorr(p_falcon,experimental));
fprintf(fresults,'%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n','Time',0,timing.gene_exp,timing.standard_fba,timing.standard_fba_best,timing.gimme,timing.shlomi,timing.fix,timing.falcon);

fclose(fresults);