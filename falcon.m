function [v_sol, corrval, nvar, v_all] = ...
    falcon(m, r, r_sd, r_group, rc, minFit, EXPCON, FDEBUG)

TESTING = false;

% If true, allow S*v >= 0 instead of S*v == 0.
MASSPROD = false;


%INPUT
%
%
%OPTIONAL INPUTS
% FDEBUG    prints or writes to files additional information
%           Important variables to always change when modifying
%           code to make sure FDEBUG is realistic:
%           NColLab <- update whenever s2 is extended
%           NRowLab <- update whenever s1 is extended
% 
%OUTPUT
%
%
%kieran: 26 apr 12
%Brandon
% 
% 
% Brandon Barker 01/15/2013  Based on Kieran Smallbone's script from:
%                            http://www.biomedcentral.com/1752-0509/6/73
%
% Improvements since first pub:
% 1) LFP to allow for an optimal scaling parameter (use initial scaling and then remove).
% 2) Parallel FVA
% 3) min disjuncti!on algorithm used for improved expression estimate and speed,
%    and better rule compatability 
% 4) Only do FVA for reactions with available gene data
% 5) Adjusted to use dual-simplex non-concurrent solver in Gurobi, since
%    problem has more rows than cols (about 2x with recon 2) and FVA is parallel.
%    This should also be the most memory-efficient algorithm in Gurobi.
% 6) Warm start FVA 
% 7) Removed repeated calls to size, which can add up to more than a second

% TODO:
%
%  Move a lot of stuff out of the while loop: this should be possible now that
%  we are using an irreversible model.
%
%  Print the reaction names/compartments that are fixed to be irrev at each
%  iteration. This may shed some light on why things sometimes go very slow.
%
%  Estimate nxmax again for spalloc
%
%  Move spalloc allocation out of while loop
%
%  Consider how to deal with massively duplicated reations: for rxns
%  with just a few duplicated gene rules, it should not be worried about much.
%  Something simple would be best. If a set of rxns has the same gene rule,
%  can we minimize the distance between the sum of absolute fluxes of the rxn
%  set and the expression value?
%  Should be possible
%  But this different for linear vs parallel "pathways". The relaxing constraint is
%  is the stringent constraint for the other, and vice vera. This could possibly be handled
%  by a trick: if the rxns have the same name, they are linear, if not, they are parallel.
%
% Remove as many find operations as possible and use logical indexing for performance:
% http://www.mathworks.com/help/matlab/math/matrix-indexing.html?s_tid=doc_12b
%
% Add fixed nvar option, e.g. for mutation simulation
%
% Use new Gurobi Tuning tool once the rest is fixed
% Use several datasets (control, couple of cancers).
%
% Determine a good upper bound for number of non-zero N values.
%
% Non-unique names for gene_to_scale
%
% Add separate module to directly call Gurobi
% 
% Add in sum of abs fluxes equals constant option. Initializing scaling
% with normalization to a particular flux may create bias, particularly
% where a long linear pathway is involved.
%
% Covariance included in standard deviation calculation (part of minDisj)
%

%
% We assume m is an irreverisble model, and when m.rev(i) == 1, 
% then m.rxn(i+1) is the reverse rxn of m.rxn(i).
%

t_falcon = tic;

if ~exist('FDEBUG', 'var')
    FDEBUG = false;
end
if ~exist('EXPCON', 'var')
    EXPCON = true;
end

nrxns = length(m.rxns);
nmets = length(m.mets);
notnan_r = ~isnan(r);
% flux_sum is used to ensure that in the LFP, we don't obtain
% a zero flux, as this is unhelpful. To estimate it, we take the
% smallest non-zero flux bound and multiply that by the number of
% irreversible reactions with expression data. In principle,
% the value shouldn't matter (as long as it is strictly positive)
% since fluxes are scaled in the LFP transform, but this is probably
% helpful for numeric stability. However, too large values will be 
% infeasible or may introduce spurious fluxes just for the sake
% of fulfulling this requirement. Too small values eventually
% cause convergence problems at later iterations: I am not sure
% if this is a numerical problem or something else. For now
% this seems to be a good approximation:

minUB = min(m.ub(m.ub > 0));
flux_sum = sum(~m.rev & notnan_r)*minUB;
%flux_sum = min(m.ub(m.ub > 0)) / 2
if flux_sum == 0
    flux_sum = mean(m.ub(m.ub > 0))/2
end
if FDEBUG
    flux_sum
end


% Typically required to be >= 0, we can require strict positivity
% due to having non-affine in our LFP.
ZMIN = 0.0001;

rgrp_notnan = r_group(notnan_r);
vbasN = [];
cbasN = [];
vbasS = [];
cbasS = [];
nSnz = sum(sum(abs(sign(m.S))));
ngroups = union(r_group, 1);
v_all = [];

ecrxns = find(any(m.rxnGeneMat, 2));
r_sum = sum(r(notnan_r));
r_pri_max = max(r);
r             = flux_sum * r / r_sum;
r_sd          = flux_sum * r_sd / r_sum;
if FDEBUG
    disp('Sum new r:');
    disp(sum(r(notnan_r)));
    disp('Max r (scaled r) is:');
    disp([r_pri_max max(r)]);
end

if TESTING
    expZtol = 2*flux_sum/nrxns
end
if FDEBUG
    r_med = median(r)
    r_min = min(r(r>0))
end

%See if this helps to get more irreversible reactions
%solFBA = optimizeCbModel(m, 'max');
%fOpt = solFBA.f;
fOpt = 0;
v_orig = zeros(nrxns, 1);
m.lb(m.c == 1) = minFit;

rxnhasgene = (sum(m.rxnGeneMat')~=0);

%rev = false(size(m.rxns));
corrval = nan;
nR_old = 0;
v_sol = zeros(size(m.rxns));
cnt = 0;
conv = 0;
nvar = nan;
fUpdate = 0;
rGrpsUsed = 0;
while sum(~m.rev) > nR_old
    cnt = cnt + 1;
    nR_old = sum(~m.rev); 
    % nnnanrev = sum((notnan_r) & m.rev) / 2;
    nnan_irr = r_group(notnan_r & ~m.rev);
    nnnan_irr = length(intersect(nnan_irr, nnan_irr));
    % Used for an alternative problem that incorporates
    % reversible reactions.
    %nnan_all = r_group(notnan_r);
    %nnnan_all = length(intersect(nnan_all, nnan_all));
    r_group_cons = zeros(1, nrxns);
    if FDEBUG
        NColLab = m.rxns;
        NRowLab = m.mets;
    end
    % 1. fit to data

    %Preallocate matrix and vectors:
    %  rows:    nmets + LFPunit + model lb/ubs + flux_sum + exp residuals + f_pre,
    %  cols:    nrxns + n + z + exp residual vars
    N = spalloc(nmets + 1       + 2*nrxns      + 1        + 2*nnnan_irr   + 1, ...
                nrxns + 1 + 1 + nnnan_irr            , floor(2.3*nSnz));
    if FDEBUG
        sz_N = size(N)
    end
    sz_N = size(N);

    N(1:nmets, 1:nrxns) = sparse(m.S);
    L = m.lb;
    U = m.ub;
    f = zeros(size(m.rxns))';
    b = zeros(size(m.mets));
    csense = '';
    if MASSPROD
        csense(1:length(b)) = 'G';
    else
        csense(1:length(b)) = 'E';
    end

    s1 = nmets; s2 = nrxns;
    %Add a column for the normalization variable
    %and the linear fractional variable z (see B&V 4.32).
    N(s1, s2 + 1) = 0; %n
    N(s1, s2 + 2) = 0; %z
    if FDEBUG
        NColLab{s2 + 1} = 'n';
        NColLab{s2 + 2} = 'z';
    end
    L(s2 + 1) = 0;
    U(s2 + 1) = inf;
    L(s2 + 2) = ZMIN;
    U(s2 + 2) = inf;
    f(s2 + 1) = 0;
    f(s2 + 2) = 0;
    s2 = s2 + 2;

    %Add unitary constraint for denominator (see B&V 4.32).
    csense(s1 + 1) = 'E';
    N(s1 + 1, nrxns + 1) = 1;
    if FDEBUG
        NRowLab{s1 + 1} = 'LFP unitary';
    end
    b(s1 + 1) = 1;
    s1 = s1 + 1;

    %Add in transformed L,U constrains (B&V 4.32).
    %consider adding conditionals here for U or L == 0.
    %since v >= 0 just gives -v =< 0 under CC transform
    %When using scaled expression to constrain vmax, use that
    %instead - this should be exclusive from medium-based
    %bounds. If not, may need to reconsider the scheme.
    for k = 1:nrxns
        f(k) = -rc; %regularization constant 
        b(s1 + 1) = 0; 
        b(s1 + 2) = 0;
        if ~EXPCON || isnan(r(k)) || U(k) == 0 % use default constraint
            N(s1 + 1, k) = 1; 
            N(s1 + 2, k) = -1; 
            N(s1 + 1, nrxns + 2) = -U(k);
            N(s1 + 2, nrxns + 2) = L(k);
        else % use expression constraint
            N(s1 + 1, k) = 1; 
            N(s1 + 2, k) = 1; 
            N(s1 + 1, nrxns + 1) = -r(k);
            N(s1 + 2, nrxns + 1) = -r(k);
        end
        L(k) = 0;
        U(k) = inf;
        if m.ub(k) == 0
            U(k) = 0;
        end
        csense(s1 + 1) = 'L';
        csense(s1 + 2) = 'L'; 
        if FDEBUG
            NRowLab{s1 + 1} = [m.rxns{k} ':U'];
            NRowLab{s1 + 2} = [m.rxns{k} ':L'];
        end
        s1 = s1 + 2;
    end 

    %Require the sum of fluxes to be above a threshold
    for k = 1:length(ecrxns)
        N(s1+1, ecrxns(k)) = -1;
    end
    flux_sum = sum(~m.rev & notnan_r)*minUB;
    N(s1 + 1, nrxns + 2) = flux_sum;
    %b(s1 + 1) = flux_sum;
    b(s1 + 1) = 0; 
    csense(s1 + 1) = 'L';
    if FDEBUG
        NRowLab{s1 + 1} = 'FlxSum';
    end
    s1 = s1 + 1;
   
    % Add a constraint on the objective value, used
    % below in: while k < nrxns
    objPriorRow = s1 + 1;
    csense(objPriorRow) = 'L';
    b(objPriorRow) = 0;
    if FDEBUG
        NRowLab{s1 + 1} = 'ObjPrior';
    end
    s1 = s1 + 1;

    k = 0;
    r_group_visited(1:nrxns) = false;
    first_r_group_visited = -1;
    rGrpsPrev = rGrpsUsed;
    rGrpsUsed = 0;
    while k < nrxns
        k = k + 1;
        d = r(k);
        s = r_sd(k);
        %objDenom = max(d, expZtol) * s;
        objDenom = s;
        cons1 = 0;
        if ~r_group_visited(r_group(k))
            cons1 = s1 + 1;
            r_group_cons(r_group(k)) = cons1;
        else
            cons1 = r_group_cons(r_group(k));
        end
        if ~m.rev(k) && ~isnan(d) && s > 0 %(s > 0 should always be true anyway)
            if first_r_group_visited <= 0
                first_r_group_visited = r_group(k);
            end
            if ~r_group_visited(r_group(k))
                r_group_visited(r_group(k)) = true;
                s1 = s1 + 2;
                rGrpsUsed = rGrpsUsed + 1; 
                if r_group(k) ~= first_r_group_visited
                    s2 = s2 + 1;
                end
            end
            %disp([k size(N,2) numel(NColLab) size(N,1) numel(b)])
            %First abs constaint:
            N(cons1, nrxns + 1) = -d;  %This is the normalization variable
            N(cons1, k) = 1;  
            N(cons1, s2 + 1) = -1;     %delta variable
            b(cons1) = 0;
            %Second abs constaint:
            N(cons1 + 1, nrxns + 1) = d; %This is the normalization variable
            N(cons1 + 1, k) = -1;  
            N(cons1 + 1, s2 + 1) = -1;  %delta variable
            b(cons1 + 1) = 0;
                    f(s2 + 1) = - 1 / objDenom;
                    L(s2 + 1) = 0;      % this can be left as 0 in the CC transform; 
                    U(s2 + 1) = inf;    % it is just the same has having -delta <= 0
            csense(cons1)   = 'L';
            csense(cons1 + 1) = 'L';
            fUpdate = fUpdate - abs(v_orig(k) - d)/objDenom;
            if FDEBUG
                NRowLab{cons1} = ['RG_' num2str(r_group(k))];
                NRowLab{cons1 + 1} = ['RG_' num2str(r_group(k))];
                NColLab{s2 + 1} = ['t_' num2str(r_group(k))];
            end    
            % Require the objective value to remain stable.
            if TESTING && cnt > 1 && abs(corrval) > 0
                N(objPriorRow, s2 + 1) = 1/objDenom;
            end

        end %end of if not nan
    end %end while k < nrxns
    
    if TESTING && cnt > 1 && abs(corrval) > 0
        N(objPriorRow, nrxns + 2) = rGrpsUsed*(corrval / rGrpsPrev);
    end

if ~all(sz_N == size(N)) || sz_N(1) ~= numel(b) || sz_N(2) ~= numel(L)
    disp('WARNING: mismatch in estimated and actual dimension detected!!!');
end

if 1
    if FDEBUG
        disp(['Not Reversible: ' num2str(sum(~m.rev))]);
        [v, fOpt, conv, vbasN, cbasN] = easyLP(f, N, b, L, U, ...
                                                  csense, vbasN, cbasN, ...
                                                  FDEBUG, NRowLab, NColLab, cnt);
    else
        [v, fOpt, conv, vbasN, cbasN] = easyLP(f, N, b, L, U, ...
                                                  csense, vbasN, cbasN);
    end
end

% This seems to do more poorly because of how the score can be defined
% by how well a single reaction matches - certainly it does bad in 
% the test model, but it needs to be tested on recon 2 to
% see if it is worth pursuing further. We also have to be careful not
% to favor longer pathways where a shorter one will do - regularization
% should help this to some extent, but again the numerics are tricky.

if 0
    revRxns = find(m.rev);
    if numel(revRxns) > 0
        firstRevRxn = revRxns(1);
        FLsave = L(firstRevRxn);
        FUsave = U(firstRevRxn);
        BLsave = L(firstRevRxn + 1);
        BUsave = U(firstRevRxn + 1);
        m.rev(firstRevRxn) = 0;
        m.rev(firstRevRxn + 1) = 0;
        % Do forward = 0 first:
        m.lb(firstRevRxn) = 0;
        m.ub(firstRevRxn) = 0;
    end
    if FDEBUG
        disp(['Not Reversible: ' num2str(sum(~m.rev))]);
        [v_b, fOpt_b, conv_b, vbasN_b, cbasN_b] = easyLP(f, N, b, L, U, ...
                                                  csense, vbasN, cbasN, ...
                                                  FDEBUG, NRowLab, NColLab, cnt);
    else
        [v_b, fOpt_b, conv_b, vbasN_b, cbasN_b] = easyLP(f, N, b, L, U, ...
                                                  csense, vbasN, cbasN);
    end

    % Do backward = 0:
    if numel(revRxns) > 0
        L(firstRevRxn) = FLsave;
        U(firstRevRxn) = FUsave;
        L(firstRevRxn + 1) = 0;
        U(firstRevRxn + 1) = 0;
    end
    if FDEBUG
        disp(['Not Reversible: ' num2str(sum(~m.rev))]);
        [v_f, fOpt_f, conv_f, vbasN_f, cbasN_f] = easyLP(f, N, b, L, U, ...
                                                  csense, vbasN, cbasN, ...
                                                  FDEBUG, NRowLab, NColLab, cnt);
    else
        [v_f, fOpt_f, conv_f, vbasN_f, cbasN_f] = easyLP(f, N, b, L, U, ...
                                                  csense, vbasN, cbasN);
    end
    L(firstRevRxn + 1) = BLsave;
    U(firstRevRxn + 1) = BUsave;
    irrev_nz_f = sum(v_f(find(~m.rev)) ~= 0)
    irrev_nz_b = sum(v_b(find(~m.rev)) ~= 0)
    nNZE_f = countNonZeroEq(v_f, m.rev, nrxns)
    nNZE_b = countNonZeroEq(v_b, m.rev, nrxns)
    if (fOpt_f / irrev_nz_f > fOpt_b / irrev_nz_b)
        [v, fOpt, conv, vbasN, cbasN] = deal(v_f, fOpt_f, conv_f, vbasN_f, cbasN_f);
        m.lb(firstRevRxn + 1) = 0;
        m.ub(firstRevRxn + 1) = 0;        
    else
        [v, fOpt, conv, vbasN, cbasN] = deal(v_b, fOpt_b, conv_b, vbasN_b, cbasN_b);
        m.lb(firstRevRxn) = 0;
        m.ub(firstRevRxn) = 0;
    end
end % of if 1/0

   
    if FDEBUG
        disp('fOpt, n, z:');
        disp([fOpt v(nrxns + 1) v(nrxns + 1)]);
    end
    if conv
        v_orig = v;
        if v(nrxns + 2) ~= 0
            v_orig = v / v(nrxns + 2); %Transform to original
            corrval = fOpt / v(nrxns + 2);
        end
        v_sol = v_orig(1:nrxns);
        v_all = [v_all v_sol];
        nvar = v_orig(nrxns + 1);
        [m.lb m.ub m.rev] = setRxnDirection(v(1:nrxns), m.lb, m.ub, m.rev, nrxns, cnt, m);
        if FDEBUG
            disp('New nvar, zvar is:');
            disp([nvar v(nrxns + 2)]);
            disp('First 15 fluxes:')
            disp(v_sol(1:15)');
            disp('Num Irrev, Previous Num Irrev:')
            disp([sum(~m.rev) nR_old]);
        end    
    end
end % end of while sum(~m.rev) > nR_old

falconTime = toc(t_falcon);
if conv
    disp(['FALCON converged in ' num2str(falconTime) ' seconds and ' ...
           num2str(cnt) ' iterations.']);
else
    disp(['FALCON did NOT converge in ' num2str(falconTime) ...
          ' seconds and ' num2str(cnt) ' iterations.']);
end
end % of falcon

function [v, fOpt, conv, svbas, scbas] = easyLP(f, a, b, vlb, vub, csense, ...
                                                vbas, cbas, FDEBUG,        ...
                                                rowLabels, colLabels, cnt)
%
%easyLP
%
% solves the linear programming problem: 
%   max f'x subject to 
%   a x = b
%   vlb <= x <= vub. 
%
% Usage: [v,fOpt,conv] = easyLP(f,a,b,vlb,vub)
%
%   f           objective coefficient vector
%   a           LHS matrix
%   b           RHS vector
%   vlb         lower bound
%   vub         upper bound
%
%   v           solution vector
%   fOpt        objective value
%   conv        convergence of algorithm [0/1]
%
% the function is a wrapper for on the "solveCobraLP" script provided with
% the COBRA (COnstraint-Based Reconstruction and Analysis) toolbox 
% http://opencobra.s.net/
%
%kieran, 20 april 2010

if ~exist('FDEBUG', 'var')
    FDEBUG = 0;
end

% matlab can crash if inputs nan
if any(isnan(f)) || any(any(isnan(a))) || any(isnan(b))...
        || any(isnan(vlb)) || any(isnan(vub)) || any(isnan(csense)) 
    error('nan inputs not allowed');
end


% initialize
v = zeros(size(vlb));
v = v(:);
f = full(f(:));
vlb = vlb(:);
vub = vub(:);

% Do this after optimization below as well.
if exist('rowLabels', 'var')
    printFalconProblem(rowLabels, colLabels, cnt, a, b, vlb, vub, f, ...
                       csense, 0);
end

% remove any tight contstraints as some solvers require volume > 0
j1 = (vlb ~= vub); 
j2 = (vlb == vub);

v(j2) = vlb(j2);
b = b(:) - a*v;
%b(isnan(b)) = 0; % likely a subtraction of nans
a(:, j2) = []; 
vlb(j2) = []; 
vub(j2) = [];
f0 = f;
f(j2) = [];
fOpt = nan;

if any(isnan(b))
    error('nan inputs not allowed: something went wrong');
end

params.method = 1;
if nargin > 6 && length(vbas) > 0
    %params.vbasis = zeros(1, length(vlb));
    %params.vbasis(1:length(vbas)) = vbas;
    params.vbasis = vbas;
    params.cbasis = zeros(1, length(b));
    params.cbasis(1:length(cbas)) = cbas;
end

if FDEBUG
    t_easy = tic; 
end
solution = solveCobraLP(...
    struct('A', a, 'b', b, 'c', f, 'lb', vlb, 'ub', vub, ...
    'osense',-1,'csense',csense) , ...
    'GurobiParams',params);
    %'printLevel',1);
if FDEBUG
    toc(t_easy)
end

% define outputs
conv = solution.stat == 1;
svbas = []; %svbas = solution.basis;
scbas = []; %scbas = solution.cbasis;

if conv
    v0 = solution.full;
    v(j1) = v0;
    fOpt = f0' * v;
    if FDEBUG
        disp(['Convergent optimum is: ' num2str(solution.obj)]);
        if isnan(fOpt)
            disp('Converged, but fOpt still nan!');
        end
        % Print again in case optimization succeeded and we can print v
        colLabels(j2) = [];
        if exist('rowLabels', 'var')
            printFalconProblem(rowLabels, colLabels, cnt + 0.1, a, b, vlb, vub, f, ...
                               csense, v(j1));
        end
    end
end
end % of easyLP



function [iLB iUB isRev] = setRxnDirection(vI, iLB, iUB, isRev, nrxns, cnt, m)
% Compute LB/UB for irrev AND rev model, as well as 
% setting the corresponding rev vectors.
%
% Scratch that, just deal with irreversible models
%
% For any forward rev reaction in an irrev model, we assume that
% its corresponding backwards reaction immediately follows it.
%
%INPUTS
% vI        Irreversible flux distribution
% rev2Irrev      Vector mapping irreversible fluxes to reversible fluxes 
%               (Generated by convertToIrreversible)

rthresh = 0.5;

%nr_rxns = length(rev2irrev);
tol = 1e-9;

k = 0;
while k < (nrxns-1) 
    k = k + 1;
    if isRev(k)
        vSum = vI(k) + vI(k+1);
        if vI(k) / vSum > rthresh
            %Forward reaction
            %rLB(i) = 0;
            iLB(k + 1) = 0;
            iUB(k + 1) = 0;

            %rrev(i) = 0;
            isRev(k) = 0;
            isRev(k + 1) = 0;
        elseif vI(k + 1) / vSum > rthresh
            %Backward reaction
            %rUB(i) = 0;
            iLB(k) = 0;
            iUB(k) = 0;

            %rrev(i) = 0;
            isRev(k) = 0;
            isRev(k + 1) = 0;      
        end
        k = k + 1;
    end
end
end % of setRxnDirection

function [iLB iUB isRev] = setFBRxnDirection(vI, iLB, iUB, isRev, nrxns, cnt, m)
% Compute LB/UB for irrev AND rev model, as well as 
% setting the corresponding rev vectors.
%
% Scratch that, just deal with irreversible models
%
% For any forward rev reaction in an irrev model, we assume that
% its corresponding backwards reaction immediately follows it.
%
%INPUTS
% vI        Irreversible flux distribution
% rev2Irrev      Vector mapping irreversible fluxes to reversible fluxes 
%               (Generated by convertToIrreversible)

rthresh = 0.5;

%nr_rxns = length(rev2irrev);
tol = 1e-8;

k = 0;
while k < (nrxns-1) 
    k = k + 1;
    if isRev(k)
        vSum = vI(k) + vI(k+1);
        if (vI(k) >= tol) && (vI(k) - vI(k+1) <= tol)
            if FDEBUG
                disp('set rxns to zero:');
                disp(m.rxns([k k+1]));
                disp(vI([k k+1]));
                disp('END set rxns to zero');
            end
            iLB(k + 1) = 0;
            iUB(k + 1) = 0;
            isRev(k) = 0;
            isRev(k + 1) = 0;
            iLB(k) = 0;
            iUB(k) = 0;
            isRev(k) = 0;
            isRev(k + 1) = 0;      
        end
        k = k + 1;
    end
end
end % of setFBRxnDirection


function nNZE = countNonZeroEq(vI, isRev, nrxns)
tol = 1e-8;
k = 0;
nNZE = 0;
while k < nrxns 
    k = k + 1;
    if isRev(k)
        if (vI(k) >= tol) && (vI(k) - vI(k+1) <= tol)
            nNZE = nNZE + 1;
        end
        k = k + 1;
    end
end
end % of countNonZeroEq