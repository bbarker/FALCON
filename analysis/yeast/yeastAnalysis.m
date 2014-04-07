function [reaction_name, experimental, p_gene_exp, p_standard_fba,        ...
    p_standard_fba_best, p_gimme, p_shlomi, p_fix, p_falcon,              ...
    s_fix, s_falcon, timing] =                                            ...
    yeastAnalysis(model, genedata_filename, experimental_fluxes_filename, ...
    gene_to_scale, flux_to_scale, methodList, nReps)

% kieran: 26 apr 12

% Brandon Barker Jan 2013 - Jan 1014

allMethods = {'FALCON', 'eMoMA', 'GIMME', 'Shlomi', 'fitFBA'};

if length(methodList) < 1
    methodList = allMethods;
end

%Whether to use new complexation method in Lee method
useMinDisj = true;
expCon = false;
minFit = 0.0;
regC = 0;

nrxns = length(model.rxns);

[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

% load transcript data
genedata	= importdata(genedata_filename);
genenames	= genedata.textdata(:,1);
genenames(1)= [];
gene_exp	= genedata.data(:,1);
gene_exp_sd	= genedata.data(:,2);

% map gene weighting to reaction weighting

expTime = tic;
[rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
    computeMinDisj(modelIrrev, genedata_filename);
minDisjTime = toc(expTime);
minDisjTime = 0; % See comment:
%
% For Yeast, this really doesn't take much time for either 
% method, so ignore expression time. For human rules, this can
% matter.

% FALCON
if find(strcmp(methodList, 'FALCON'))
    disp('Running FALCON ...');
    %[rxn_exp, rxn_exp_sd, rxn_missing_gene] = geneToRxn(model, genedata_filename);

    [rxn_exp_irr, rxn_exp_sd_irr] = geneToReaction(modelIrrev, genenames, ...
        gene_exp, gene_exp_sd);
    % sds 0 -> small
    rxn_exp_sd_irr(rxn_exp_sd_irr == 0) = min(rxn_exp_sd_irr(rxn_exp_sd_irr>0))/2;
    %First compare the number of nans and zeros in both data typs.
    %gte_nan = sum(isnan(rxn_exp_irr))
    %mdj_nan = sum(isnan(rxn_exp_md))
    %gte_0 = sum(rxn_exp_irr == 0)
    %mdj_0 = sum(rxn_exp_md == 0)
    rxnRuleCell = cell(length(model.rxns) + 1, 3);
    rxnRuleCell{1,1} = 'Rule';
    rxnRuleCell{1,2} = 'geneToRxn';
    rxnRuleCell{1,3} = 'minDisj';
    for i = 1:length(model.rxns)
        rxnRuleCell{i + 1, 1} = modelIrrev.grRules{i};
        rxnRuleCell{i + 1, 2} = num2str(rxn_exp_irr(i));
        rxnRuleCell{i + 1, 3} = num2str(rxn_exp_md(i));
    end
    cell2csv(['ComplexExpressionCompare.csv'], rxnRuleCell, ',');

    tic;
    %Need to separate transcript data loading
    if useMinDisj % need to modify to use old expression values
        [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s] =         ...
            falconMulti(modelIrrev, nReps, rxn_exp_md,        ...
            rxn_exp_sd_md, rxn_rule_group, 'rc', regC,        ...
            'minFit', minFit, 'EXPCON', expCon);
    else
        [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s] =          ...
            falconMulti(modelIrrev, nReps,  rxn_exp_irr,       ...
            rxn_exp_sd_irr, rxn_rule_group, 'rc', regC,        ...
            'minFit', minFit, 'EXPCON', expCon);
    end
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);
    timing.falcon = toc/nReps + minDisjTime;
    v_falcon_s = convertIrrevFluxDistribution(v_falconIrr_s, matchRev);
    save([genedata_filename '_falcon_flux.mat'], 'v_falcon');
else
    v_falcon = zeros(length(model.lb), 1);
    v_falcon_s = v_falcon;
    timing.falcon = 0 + minDisjTime;
end


% Standard FBA
% Always run FBA ... why not.
tic;
solution        = optimizeCbModel(model,[],'one');
v_standard_fba	= solution.x;
timing.standard_fba = toc;
fOpt            = solution.f;


if find(strcmp(methodList, 'GIMME'))
    m = model;
    if ~exist('rxn_exp', 'var')
        [rxn_exp, rxn_exp_sd] = geneToReaction(m, genenames, ...
            gene_exp, gene_exp_sd);
    end
    disp('Running GIMME ...');
    % "required metabolic functionalities" (=growth) set to 90% of maximum
    m.lb(m.c == 1)	= 0.9*fOpt;
    % gimme
    tic;
    v_gimme         = gimme(m, rxn_exp);
    timing.gimme = toc;
    save([genedata_filename '_gimme_flux.mat'], 'v_gimme');
else
    v_gimme = zeros(length(model.lb), 1);
    timing.gimme = 0;
end

% shlomi
if find(strcmp(methodList, 'Shlomi'))
    if ~exist('rxn_exp', 'var')
        [rxn_exp, rxn_exp_sd] = geneToReaction(model, genenames, ...
            gene_exp, gene_exp_sd);
    end
    disp('Running Shlomi ...');
    tic;
    v_shlomi    	= shlomi(model, rxn_exp);
    timing.shlomi = toc;
    save([genedata_filename '_shlomi_flux.mat'], 'v_shlomi');
else
    v_shlomi = zeros(length(model.lb), 1);
    timing.shlomi = 0;
end

if find(strcmp(methodList, 'eMoMA'))
    m = model;
    if ~exist('rxn_exp', 'var')
        [rxn_exp, rxn_exp_sd] = geneToReaction(m, genenames, ...
            gene_exp, gene_exp_sd);
    end
    % scale by uptake reaction
    uptake              = find(strcmp(gene_to_scale,m.rxnNames));
    rxn_exp_sd          = rxn_exp_sd/rxn_exp(uptake);
    rxn_exp             = rxn_exp/rxn_exp(uptake);

    m.lb(uptake)	= 1;
    m.ub(uptake)	= 1;
    % Gene expression constraint FBA
    disp('Running eMoMA (original) ...');
    tic;
    %if useMinDisj
    %    v_gene_exp = dataToFlux(m, rxn_exp_md, rxn_exp_sd_md);
    %else
        v_gene_exp = dataToFlux(m, rxn_exp, rxn_exp_sd);
    %end
    timing.gene_exp=toc;
    save([genedata_filename '_gene_exp_flux.mat'], 'v_gene_exp');
    % fixed expression method
    disp('Running eMoMA (exp fix) ...');
    tic;
    %if useMinDisj
    %    v_fix = dataToFluxFix(m, rxn_exp_md, rxn_exp_sd_md);
    %    [v_fix, ~, v_fix_s, ~] = dataToFluxFixMulti( ...
    %        m, rxn_exp_md, rxn_exp_sd_md, nReps)
    %else
        % v_fix = dataToFluxFix(m, rxn_exp, rxn_exp_sd);
        [v_fix, ~, v_fix_s, ~] = dataToFluxFixMulti( ...
            m, rxn_exp, rxn_exp_sd, nReps)
    %end
    timing.fix = toc/nReps;
    save([genedata_filename '_fix_flux.mat'], 'v_fix');
else
    v_gene_exp = zeros(length(model.lb), 1);
    timing.gene_exp = 0;
    v_fix = zeros(length(model.lb), 1);
    v_fix_s = v_fix;
    timing.fix = 0;
end

% compare
experimental_fluxes = importdata(experimental_fluxes_filename);

reaction_name   = experimental_fluxes.textdata;
experimental	= zeros(size(experimental_fluxes.textdata,1),1);
p_gene_exp      = zeros(size(experimental_fluxes.textdata,1),1);
p_standard_fba  = p_gene_exp;
p_gimme         = p_gene_exp;
p_shlomi     	= p_gene_exp;
p_fix     	= p_gene_exp;
p_falcon     	= p_gene_exp;
s_fix           = p_gene_exp;
s_falcon        = p_gene_exp;

flux = strcmp(flux_to_scale,reaction_name);
flux = experimental_fluxes.data(flux,1);

for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(reaction_name{k}, model.rxnNames));
    experimental(k)     = experimental_fluxes.data(k,1);
    if k == 1
         experimental(k) = -experimental(k);
    end
    p_gene_exp(k)       = flux*v_gene_exp(j);
    p_standard_fba(k)	= flux*v_standard_fba(j);
    p_gimme(k)          = flux*v_gimme(j);
    p_shlomi(k)     	= flux*v_shlomi(j);
    p_fix(k)     	= flux*v_fix(j);
    p_falcon(k)     	= flux*v_falcon(j);
    s_fix(k)     	= flux*v_fix_s(j);
    s_falcon(k)     	= flux*v_falcon_s(j);
    % reaction_name{k}
    % v_falcon(j)
end

% remove small entries
p_gene_exp(abs(p_gene_exp)<1e-6)            = 0;
p_standard_fba(abs(p_standard_fba)<1e-6)	= 0;
p_gimme(abs(p_gimme)<1e-6)                  = 0;
p_shlomi(abs(p_shlomi)<1e-6)                = 0;
p_fix(abs(p_fix)<1e-6)                = 0;
p_falcon(abs(p_falcon)<1e-6)                = 0;

% find best fit from standard FBA solution
% ... overkill for this problem, but reuses existing method
if find(strcmp(methodList, 'fitFBA'))
    tic;
    model.lb(model.c == 1) = fOpt;
    data = nan(size(model.rxns));
    uptake              = find(strcmp(gene_to_scale,model.rxnNames));
    data(uptake) = 1;
    for k = 1:size(experimental_fluxes.textdata,1)
        j = find(strcmp(experimental_fluxes.textdata{k,1},model.rxnNames));
        data(j) = experimental_fluxes.data(k,1)/flux; %#ok<FNDSB>
    end
    v_standard_fba_best = dataToFlux(model,data,data);
    timing.standard_fba_best=toc;
else
    v_standard_fba_best = zeros(length(model.lb), 1);
    timing.standard_fba_best = 0;
end
p_standard_fba_best	= zeros(size(experimental_fluxes.textdata,1),1);
for k = 1:size(experimental_fluxes.textdata,1)
    j = find(strcmp(experimental_fluxes.textdata{k,1},model.rxnNames));
    p_standard_fba_best(k) = flux*v_standard_fba_best(j); %#ok<FNDSB>
end
p_standard_fba_best(abs(p_standard_fba_best)<1e-6)	= 0;


function v_sol = dataToFlux(m,r,r_sd)

% kieran: 21 sep 11

rev = false(size(m.rxns));

nR_old = 0;

% m.lb(~ismember(m.lb,[-inf,0,inf,-1000,1000])) = -inf;
% m.ub(~ismember(m.ub,[-inf,0,inf,-1000,1000])) = inf;

v_sol = zeros(size(m.rxns));

while sum(~m.rev) > nR_old
    
    nR_old = sum(~m.rev);
    
    % 1. fit to data
    
    N = m.S;    
    L = m.lb;
    U = m.ub;
    f = zeros(size(m.rxns))';
    b = zeros(size(m.mets));
    
    for k = 1:length(m.rxns)
        d = r(k);
        s = r_sd(k);
        
        if ~m.rev(k) && ~isnan(d) && s>0
            [s1,s2] = size(N);
            N(s1+1,k) = 1; N(s1+1,s2+1) = -1; N(s1+1,s2+2) = 1;
            L(s2+1) = 0; L(s2+2) = 0;
            U(s2+1) = inf; U(s2+2) = inf;
            b(s1+1) = d;
            f(s2+1) = - 1/s;
            f(s2+2) = - 1/s;
        end
    end
    
    [v,fOpt,conv] = easyLPLee(f,N,b,L,U);
    
    if conv
        v_sol = v(1:length(m.rxns));
        for k = 1:length(m.rxns)
            if rev(k), v_sol(k) = - v_sol(k); end
        end
        
        % 2. run FVA
        
        N = [N; f]; %#ok<AGROW>
        b = [b(:); fOpt];
        
        for k = 1:length(m.rxns)
            
            if m.rev(k)
                
                f = zeros(size(L));
                f(k) = -1;
                
                [~,fOpt,conv] = easyLPLee(f,N,b,L,U);
                
                if conv && (-fOpt >= 0) % irreversibly forward
                    
                    m.lb(k) = max(m.lb(k),0);
                    m.rev(k) = 0;
                    
                else
                    f(k) = 1;
                    [~,fOpt,conv] = easyLPLee(f,N,b,L,U);
                    
                    if conv && abs(fOpt)<=0 % irreversibly backward
                        
                        m.S(:,k) = - m.S(:,k);
                        
                        m.ub(k) = - m.ub(k);
                        m.lb(k) = - m.lb(k);
                        
                        ub = m.ub(k);
                        m.ub(k) = m.lb(k);
                        m.lb(k) = ub;
                        
                        m.lb(k) = max(m.lb(k),0);
                        m.rev(k) = 0;
                        
                        rev(k) = ~rev(k);
                        
                    end
                end
            end
        end
    end
end




function v_sol = gimme(model,gene_exp)

% cutoff at lower quartile
cutoff                  = quantile(gene_exp(~isnan(gene_exp)),0.25);

% set up big matrices
model.c = zeros(size(model.c));

for k = 1:length(gene_exp)
    if gene_exp(k) < cutoff
        c = cutoff - gene_exp(k);
        [n1,n2] = size(model.S);
        % v = v+ - v-;
        model.S(n1+1,k)    = 1; model.b(n1+1) = 0;
        model.S(n1+1,n2+1) = -1; model.lb(n2+1) = 0; model.ub(n2+1) = inf; model.c(n2+1) = -c;
        model.S(n1+1,n2+2) = 1; model.lb(n2+2) = 0; model.ub(n2+2) = inf; model.c(n2+2) = -c;
    end
end

solution        = optimizeCbModel(model,[],'one');
if numel(solution.x) < 1
    solution.x = zeros(size(model.lb));
end
v_sol = solution.x(1:length(model.rxns));

function v_sol = shlomi(model,rxn_exp)

% cutoffs at lower and upper quantiles
ExpressedRxns = model.rxns(rxn_exp >= quantile(rxn_exp(~isnan(rxn_exp)),0.75));
UnExpressedRxns = model.rxns(rxn_exp <= quantile(rxn_exp(~isnan(rxn_exp)),0.25));

% below copied from createTissueSpecificModel.m, part of the Cobra toolbox
% http://opencobra.sf.net/

RHindex = findRxnIDs(model,ExpressedRxns);
RLindex = findRxnIDs(model,UnExpressedRxns);

S = model.S;
lb = model.lb;
ub = model.ub;
eps = 1;

% Creating A matrix
A = sparse(size(S,1)+2*length(RHindex)+2*length(RLindex),size(S,2)+2*length(RHindex)+length(RLindex));
[m,n,s] = find(S);
for i = 1:length(m)
    A(m(i),n(i)) = s(i); %#ok<SPRIX>
end

for i = 1:length(RHindex)
    A(i+size(S,1),RHindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1),i+size(S,2)) = lb(RHindex(i)) - eps; %#ok<SPRIX>
    A(i+size(S,1)+length(RHindex),RHindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1)+length(RHindex),i+size(S,2)+length(RHindex)+length(RLindex)) = ub(RHindex(i)) + eps; %#ok<SPRIX>
end

for i = 1:length(RLindex)
    A(i+size(S,1)+2*length(RHindex),RLindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1)+2*length(RHindex),i+size(S,2)+length(RHindex)) = lb(RLindex(i)); %#ok<SPRIX>
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),RLindex(i)) = 1; %#ok<SPRIX>
    A(i+size(S,1)+2*length(RHindex)+length(RLindex),i+size(S,2)+length(RHindex)) = ub(RLindex(i)); %#ok<SPRIX>
end

% Creating csense
csense1(1:size(S,1)) = 'E';
csense2(1:length(RHindex)) = 'G';
csense3(1:length(RHindex)) = 'L';
csense4(1:length(RLindex)) = 'G';
csense5(1:length(RLindex)) = 'L';
csense = [csense1 csense2 csense3 csense4 csense5];

% Creating lb and ub
lb_y = zeros(2*length(RHindex)+length(RLindex),1);
ub_y = ones(2*length(RHindex)+length(RLindex),1);
lb = [lb;lb_y];
ub = [ub;ub_y];

% Creating c
c_v = zeros(size(S,2),1);
c_y = ones(2*length(RHindex)+length(RLindex),1);
c = [c_v;c_y];

% Creating b
b_s = zeros(size(S,1),1);
lb_rh = lb(RHindex);
ub_rh = ub(RHindex);
lb_rl = lb(RLindex);
ub_rl = ub(RLindex);
b = [b_s;lb_rh;ub_rh;lb_rl;ub_rl];

% Creating vartype
vartype1(1:size(S,2),1) = 'C';
vartype2(1:2*length(RHindex)+length(RLindex),1) = 'B';
vartype = [vartype1;vartype2];

MILPproblem.A = A;
MILPproblem.b = b;
MILPproblem.c = c;
MILPproblem.lb = lb;
MILPproblem.ub = ub;
MILPproblem.csense = csense;
MILPproblem.vartype = vartype;
MILPproblem.osense = -1;
MILPproblem.x0 = [];

solution = solveCobraMILP(MILPproblem);

x = solution.cont;
for i = 1:length(x)
    if abs(x(i)) < 1e-6
        x(i,1) = 0;
    end
end

v_sol = x;
