function randomExpressionSim(model, expFile, nSims, label, ...
    experimental_fluxes_filename)

%Make sure to test this with just ~16 simulations before going all out.

expCon = false;
minFit = 0.0;
regC = 0;
flux_to_scale   = 'D-glucose exchange';
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

corrTypes = {'Pearson', 'Kendall', 'Spearman'};
pCorr = @(f,y)(corr(f, y, 'type', 'Pearson'));
kCorr = @(f,y)(corr(f, y, 'type', 'Kendall'));
sCorr = @(f,y)(corr(f, y, 'type', 'Spearman'));

nCTypes = length(corrTypes);
PearsonVec     = zeros(nSims, 1);
KendallVec     = zeros(nSims, 1);
SpearmanVec    = zeros(nSims, 1);
PearsonExpVec  = zeros(nSims, 1);
KendallExpVec  = zeros(nSims, 1);
SpearmanExpVec = zeros(nSims, 1);

%get expression file name
[pathstr, expName, ext] = fileparts(expFile);

genedata         = importdata(expFile);
[ndrows, ndcols] = size(genedata.data);
if ndcols == 2
    genenames	= genedata.textdata(:,1);
    genenames(1)= [];
    gene_exp	= genedata.data(:,1);
    gene_exp_sd	= genedata.data(:,2);
elseif ndcols == 3
    %genenames	= cellfun(@num2str, ...
    %    columnVector(num2cell(genedata.data(:, 1))), 'UniformOutput', false);
    genenames    = genedata.data(:,1); 
    %genenames(1)= [];
    gene_exp	 = genedata.data(:,2);
    gene_exp_sd	 = genedata.data(:,3);
end

ngenes = length(gene_exp);

if ~exist('label', 'var')
    label = '';
end

%Compute unperturbed for reference.
[rxn_expRef, rxn_exp_sdRef, rxn_rule_groupRef] = ...
    computeMinDisj(modelIrrev, expFile);
[v_falconIrrRef, corrval_falconRef] = falcon(modelIrrev,      ...
    rxn_expRef, rxn_exp_sdRef, rxn_rule_groupRef, 'rc', regC, ...
    'minFit', minFit, 'EXPCON', expCon);
v_falconRef = convertIrrevFluxDistribution(v_falconIrrRef, matchRev);
disp(['v_fNormRef = ' num2str(norm(v_falconRef, 1))]);

% compare
if exist('experimental_fluxes_filename', 'var')
    experimental_fluxes = importdata(experimental_fluxes_filename);
    reaction_name   = experimental_fluxes.textdata;
    flux = strcmp(flux_to_scale, reaction_name);
    flux = experimental_fluxes.data(flux, 1);
    p_falcon      = zeros(size(experimental_fluxes.textdata,1),1);
    for k = 1:size(experimental_fluxes.textdata,1)
        j = find(strcmp(reaction_name{k},model.rxnNames));
        experimental(k)     = experimental_fluxes.data(k,1);
        p_falcon(k)         = flux*v_falconRef(j);
    end
    experimental = columnVector(experimental);

    % Kieran said the signs from the Lee paper flux data matched
    % the signs predicted by the method; this equates to
    % only glucose being negative:
    experimental(1) = -1*experimental(1);

    pCorrRef  = pCorr(p_falcon, experimental)
    kCorrRef  = kCorr(p_falcon, experimental)
    sCorrRef = sCorr(p_falcon, experimental)

end

startTime = strrep(strrep(num2str(clock()), ' ', ''), '.', '');

parfor i = 1:nSims
    permVec = randperm(ngenes);
    gene_exp_perm = gene_exp(permVec);
    gene_exp_sd_perm = gene_exp_sd(permVec);
    %genenames_perm = genenames(permVec);
    % Since FALCON currently only supports calling minDisj on a file,
    % this means we need to write out all the data to a randomly
    % generated file and then delete the file.
    tmpFileName = [label 'randCorr_' num2str(nSims) '_' num2str(i) ...
                   '_' expName '_' startTime '.csv'];
    gdout = genedata;
    if ndcols == 2
        % need permuted gene names too
        % gdout.textdata(2:end, 1) = genenames_perm;
        gdout.textdata(2:end, 2) = cellfun(@num2str, ...
            num2cell(gene_exp_perm), 'UniformOutput', false);
        gdout.textdata(2:end, 3) = cellfun(@num2str, ...
            num2cell(gene_exp_sd_perm), 'UniformOutput', false);
        cell2csv(tmpFileName, gdout.textdata, '\t');
    elseif ndcols == 3
        %gdout.data(:,1) = genenames_perm;
        gdout.data(:,2) = gene_exp_perm;
        gdout.data(:,3) = gene_exp_sd_perm; 
        cell2csv(tmpFileName, gdout.textdata(1,:), '\t');
        dlmwrite(tmpFileName, gdout.data, '-append', 'delimiter', ...
            '\t', 'precision', 15);
    else
        disp('Problem reading gene expression file.');
        %return; ?? what to do in parfor?
    end    
    [rxn_exp, rxn_exp_sd, rxn_rule_group] = ...
        computeMinDisj(modelIrrev, tmpFileName);
    [v_falconIrr, corrval_falcon] = falcon(modelIrrev,    ...
        rxn_exp, rxn_exp_sd, rxn_rule_group, 'rc', regC,  ...
        'minFit',  minFit, 'EXPCON', expCon);
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);
    delete(tmpFileName);
%    if exist('experimental_fluxes_filename', 'var') 
%    No way to test for existence in parfor?
        %compute experimental correlation
        p_falcon = zeros(size(experimental_fluxes.textdata, 1), 1);


        for k = 1:size(experimental_fluxes.textdata,1)
            j = find(strcmp(reaction_name{k}, model.rxnNames));
            p_falcon(k)     	= flux*abs(v_falcon(j));
        end
        % remove small entries
        p_falcon(abs(p_falcon) < 1e-6) = 0;
        PearsonExpVec(i)  = pCorr(p_falcon, experimental);
        KendallExpVec(i)  = kCorr(p_falcon, experimental);
        SpearmanExpVec(i) = sCorr(p_falcon, experimental);
%    end

    PearsonVec(i)  = pCorr(v_falcon, v_falconRef);
    KendallVec(i)  = kCorr(v_falcon, v_falconRef); 
    SpearmanVec(i) = sCorr(v_falcon, v_falconRef);
end

outFileName = [label 'randCorr_' num2str(nSims) '_' expName '.mat'];

save(outFileName, 'PearsonExpVec', 'KendallExpVec', 'SpearmanExpVec', ...
     'PearsonVec', 'KendallVec', 'SpearmanVec');
