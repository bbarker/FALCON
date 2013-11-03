function [rxn_exp,rxn_exp_sd,rxn_rule_group] = ...
    computeMinDisj(model, genedata_filename, sigma, FDEBUG)
% Requires the cell2csv package (for now; need to change to FIFOs)
% minDisj needs to be in $PATH (system path)

if ~exist('FDEBUG', 'var')
    FDEBUG = false;
end
if ~exist('sigma', 'var')
    sigma = 0;
end

ztol = 1e-4;

mdesc = strrep(model.description, ' ', '_');
rfid = num2str(randint(1,1,10e40));
rfname = [genedata_filename, '_', mdesc, rfid];
rfname = strrep(rfname,' ', '');
cell2csv(rfname, model.grRules,',',2000);
rfout = [rfname, '_out'];
nrxns = length(model.rxns);
if FDEBUG
    disp('Filebase: ');
    disp(rfname);
end
%Perturb the expression vector
if sigma > 0
    gfr_fid = num2str(randint(1,1,10e40));
    genedata_filename_pert = [genedata_filename, '_', gfr_fid]
    genedata = importdata(genedata_filename, '\t')  
    [ndrows, ndcols] = size(genedata.data);
    randVec = lognrnd(-sigma^2/2,sigma,ndrows,1);
    if ndcols == 2
        %mede = median(genedata.data(~isnan(genedata.data(:,1)),1));
        %randVec = simpleTruncatedNorm(sigma, 0, inf, ndrows, mede);
        pert_vec = genedata.data(:,1) .* randVec;
        %cell2csv on textdata after perturbing first data column
        genedata.textdata(2:end,2) = cellfun(@num2str, ...
            num2cell(pert_vec), 'UniformOutput', false);
        genedata.textdata(2:end,3) = cellfun(@num2str, ...
            num2cell(genedata.data(:,2)), 'UniformOutput', false);
        cell2csv(genedata_filename_pert, genedata.textdata, '\t', 2000);
    elseif ndcols == 3
        %mede = median(genedata.data(~isnan(genedata.data(:,2)),2))
        %randVec = simpleTruncatedNorm(sigma, 0, inf, ndrows, mede);
        %randVec = randVec(:);
        genedata.data(:,2) = genedata.data(:,2) .* randVec; 
        cell2csv(genedata_filename_pert, genedata.textdata(1,:), '\t', 2000);
        dlmwrite(genedata_filename_pert, genedata.data, '-append', 'delimiter', ...
            '\t', 'precision', 15);
    else
        disp('Problem reading gene expression file.');
        return;
    end
    genedata_filename = genedata_filename_pert;
end

[status, cmdout] = system(['minDisj ', genedata_filename, ' ', rfname, ' > ', rfout]);
if status ~= 0
    pause(0.03); %why (or) is this necessary?
    disp('try #2...');
    [status, cmdout] = system(['minDisj ', genedata_filename, ' ', rfname, ' > ', rfout]);
    if status ~= 0
        pause(3); %why (or) is this necessary?
        disp('try #3...');
        [status, cmdout] = system(['minDisj ', genedata_filename, ' ', rfname, ' > ', rfout]);
        if status ~= 0 
            disp(['minDisj failed with return code ' num2str(status)]);
            return;
        end
    end
end

if FDEBUG == 0
    delete(rfname)
end
if sigma > 0 && FDEBUG == 0
    delete(genedata_filename)
end

%create rxn_rule_group
%It may be better to allocate the appropriately sized Map first
%using set functions to get the right size.
ruleFirstIdx = containers.Map;
cidx = 0;
rxn_rule_group = zeros(nrxns, 1);
while cidx < nrxns
    cidx = cidx + 1;
    key = strtrim(model.grRules{cidx});
    if length(key) > 0
        if isKey(ruleFirstIdx,key)
            rxn_rule_group(cidx) = ruleFirstIdx(key);
        else
            rxn_rule_group(cidx) = cidx;
            ruleFirstIdx(key) = cidx;
        end
    else
        rxn_rule_group(cidx) = cidx;
    end
end

rxnData = importdata(rfout);
if FDEBUG == 0
    delete(rfout)
end
%Handle absent newline at last-line issue
if length(rxnData) == (nrxns-1)
    rxnData(nrxns,:) = nan*ones(1,2);
end

% map gene weighting to reaction weighting
%[rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);

rxn_exp = rxnData(:,1);
rxn_exp_sd = rxnData(:,2); %actually var
%We don't use this anymore:
%rxn_missing_gene = isnan(rxn_exp);
% sds 0 -> small
if FDEBUG
    disp('min max (rxn_exp then rxn_exp_sd)');
    disp([min(rxn_exp) max(rxn_exp)]);
    disp([min(rxn_exp_sd) max(rxn_exp_sd)]);
end
if max(rxn_exp_sd) > ztol
    rxn_exp_sd(rxn_exp_sd < ztol | isnan(rxn_exp_sd)) = ...
         min(rxn_exp_sd(rxn_exp_sd >= ztol))/2;
else 
    rxn_exp_sd(rxn_exp_sd == 0 | isnan(rxn_exp_sd)) = 1;
end
if FDEBUG
    disp('New min rxn_exp_sd');
    disp(min(rxn_exp_sd));
end




