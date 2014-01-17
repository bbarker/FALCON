function falconCorr = yeastResults(model, methodList, experiment, writeRes, fLabel, nReps)

% kieran: 26 apr 12
% Brandon Barker 01/19/13 

  timeInit = num2str(now());

%model originally used from paper:
%model_name      = 'yeast_5.21_MCISB.xml';
gene_to_scale   = 'glucose transport';
flux_to_scale   = 'D-glucose exchange';
%model           = readCbModel(model_name);

corrType = 'Pearson';
mycorr = @(f,y)(corr(f, y, 'type', corrType));

expFiles = {};
fluxFiles = {};
for i = 1:length(experiment)
    if experiment(i) == 75
        expFiles{end + 1} = 'genedata_75.txt';
        fluxFiles{end + 1} = 'experimental_fluxes_75.txt';
    end
    if experiment(i) == 85
        expFiles{end + 1} = 'genedata_85.txt';
        fluxFiles{end + 1} = 'experimental_fluxes_85.txt';
    end
end

if ~exist('writeRes', 'var')
    writeRes = true;
end
if ~exist('nReps', 'var')
    nReps = 7; %Reasonable default
end


falconCorr = -inf;

for i = 1:length(expFiles)
    [reaction_name, experimental, p_gene_exp, p_standard_fba,               ...
        p_standard_fba_best, p_gimme, p_shlomi, p_fix, p_falcon,            ...
        s_fix, s_falcon, timing] =                                          ...
        yeastAnalysis(model, expFiles{i}, fluxFiles{i},                     ...
            gene_to_scale, flux_to_scale, methodList, nReps);

    falconCorr = mycorr(p_falcon, experimental);

    if writeRes
        fresults = fopen([expFiles{i} '_results_' fLabel '_' timeInit '.csv'],'w');

        fprintf(fresults, ['%s' repmat('\t%s', 1, 10) '\n'],                  ...
            'Reaction', 'Experimental%', 'Expression', 'Standard FBA',       ...
	    'Fitted FBA', 'Gimme', 'iMAT', 'Expression (FVA Fixed)',         ...
            'STD: Expression (FVA Fixed)', 'falcon', 'STD: falcon');
        for k = 1:size(reaction_name,1)
            fprintf(fresults, ['%s' repmat('\t%g', 1, 10) '\n'],                        ...
                reaction_name{k,1}, experimental(k), p_gene_exp(k), p_standard_fba(k), ...
		p_standard_fba_best(k), p_gimme(k), p_shlomi(k), p_fix(k), s_fix(k),   ...
                p_falcon(k), s_falcon(k));
        end

        fprintf(fresults, ['%s' repmat('\t%g', 1, 10) '\n'], corrType, 1,              ...
            mycorr(p_gene_exp, experimental), mycorr(p_standard_fba, experimental),   ...
            mycorr(p_standard_fba_best,experimental), mycorr(p_gimme, experimental),  ...
            mycorr(p_shlomi, experimental), mycorr(p_fix, experimental), nan,       ...
	    falconCorr, nan);

        fprintf(fresults, ['%s' repmat('\t%g', 1, 10) '\n'], 'Time', 0,       ...
            timing.gene_exp, timing.standard_fba, timing.standard_fba_best,  ...
	    timing.gimme, timing.shlomi, timing.fix, nan, timing.falcon, nan);

        fclose(fresults);
    end
end % of for i = 1:length(expFiles)

