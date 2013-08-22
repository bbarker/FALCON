%kieran: 26 apr 12
%Brandon
% 
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
% Check to see if it is better to start with splitTurn = false
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
function [v_sol, corrval, nvar, v_all] = falcon(m,r,r_sd,r_group,flux_sum,rc,minFit)
%
%We assume m is an irreverisble model, and w!hen m.rev(i) == 1, then m.rxn(i+1) is the reverse rxn of m.rxn(i).
%
nrxns = length(m.rxns);
nmets = length(m.mets);
%nnnan = sum(~isnan(r));
nnnan = r_group(~isnan(r));
nnnan = length(intersect(nnnan,nnnan));
vbasN = [];
cbasN = [];
vbasS = [];
cbasS = [];
nSnz = sum(sum(abs(sign(m.S))));
ngroups = union(r_group,1);
v_all = [];

% Brandon Barker 01/15/2013
%upt_ub_save = zeros(1, length(upt_const));
%upt_lb_save = zeros(1, length(upt_const));
%for i = 1:length(gene_to_scale)
%  uptake              = find(strcmp(gene_to_scale{i},m.rxnNames));
%  upt_lb_save(i)        = m.lb(uptake)
%  m.lb(uptake)	= upt_const{i};
%  upt_ub_save(i)        = m.ub(uptake)
%  m.ub(uptake)	= upt_const{i};  
%end
%Scale expression by first uptake flux argument
%in the first iteration:
%nvar = upt_const{1};
%nvarp = nvar;
%uptake        = find(strcmp(gene_to_scale{1},m.rxnNames));
ecrxns = find(any(m.rxnGeneMat,2));
r_sum = sum(r(~isnan(r)));
r_pri_max = max(r);
r             = flux_sum * r/r_sum;
r_sd          = flux_sum * r_sd/r_sum;
disp('Sum new r:');
disp(sum(r(~isnan(r))));


disp('Max r (scaled r) is:');
disp([r_pri_max max(r)]);

%See if this helps to get more irreversible reactions
solFBA = optimizeCbModel(m,'max');
fOpt = solFBA.f;
m.lb(m.c == 1)	= minFit; %.1249

infval = inf;
% kieran: 21 sep 11

%grxns = find(sum(m.rxnGeneMat'));
rxnhasgene = (sum(m.rxnGeneMat')~=0);

%rev = false(size(m.rxns));
corrval = 0;
nR_old = 0;
v_sol = zeros(size(m.rxns));
%firstRun = true
cnt = 0;
while sum(~m.rev) > nR_old
    cnt = cnt+1;
    nR_old = sum(~m.rev); 
    nnnanrev = sum((~isnan(r)) & m.rev)/2;
    nnnan_nrev = sum((~isnan(r)) & ~m.rev);
    r_group_cons = zeros(1,nrxns);
    % 1. fit to data

    %Preallocate matrix and vectors:
    %Eventually need to test this to see how much i'm off by for 3rd argument
    %Sparcity estimate may need to be tweaked, check model:
    N = spalloc(nmets + 1 + 2*nnnan + 1, nrxns + 2 + nnnan_nrev, floor(2.1*nSnz));
    disp('size N:');
    disp(size(N));
    N(1:nmets,1:nrxns) = sparse(m.S);
    %Q = spalloc(nrxns+5*nnnanrev+2*nnnan_nrev, nrxns+5*nnnanrev+2*nnnan_nrev,nnnanrev); 
    L = m.lb;
    U = m.ub;
    disp('L init'); disp(size(L));
    f = zeros(size(m.rxns))';
    b = zeros(size(m.mets));
    csense = '';
    csense(1:length(b)) = 'E';

    s1 = nmets; s2 = nrxns;
    %if ~splitTurn
      %Add a column for the normalization variable
      %and the linear fractional variable z (see B&V 4.32).
      N(s1,s2+1) = 0; %n
      N(s1,s2+2) = 0; %z
      L(s2+1) = -infval;
      U(s2+1) = infval;
      L(s2+2) = 0;
      U(s2+2) = infval;
      f(s2+1) = 0;
      f(s2+2) = 0;
      s2 = s2+2;
 
      %Add unitary constraint for denominator (see B&V 4.32).
      csense(s1+1) = 'E';
      N(s1+1, nrxns+1) = 1;
      b(s1+1) = 1;
      s1 = s1+1;

      %Add in transformed L,U constrains (B&V 4.32).
      %consider adding conditionals here for U or L == 0.
      %since v >= 0 just gives -v =< 0 under CC transform
      for k = 1:nrxns
	f(k) = -rc; %regularization constant 
        N(s1+1, k) = 1; N(s1+1, nrxns+2) = -U(k);
        b(s1+1) = 0; csense(s1+1) = 'L';
        N(s1+2, k) = -1; N(s1+2, nrxns+2) = L(k);
        b(s1+2) = 0; csense(s1+2) = 'L'; 
        L(k) = -infval; U(k) = infval;
	s1 = s1+2;
      end 
    %Require the sum of fluxes to be above a threshold
    %It is worth noting that these are the
    %enzymatic fluxes to keep things more appropriately scaled
    for k = 1:length(ecrxns)
      N(s1+1, ecrxns(k)) = -1;
    end
    N(s1+1, nrxns+2) = flux_sum;
    b(s1+1) = 0; csense(s1+1) = 'L';
    s1 = s1+1;
 
    %end % end of if not splitTurn
    k = 0;
    while k < nrxns
        k = k+1;
        d = r(k);
        s = r_sd(k);
	cons1 = 0;
	if k == r_group(k)
	  cons1 = s1+1;
	  r_group_cons(k) = cons1;
	else
	  cons1 = r_group_cons(r_group(k));
	end
        if ~isnan(d) %&& s>0 %(assumed)
	    %First abs constaint:
	    N(cons1,nrxns+1) = -d; %This is the normalization variable
	    N(cons1,k) = 1;  
            N(cons1,s2+1) = -1; %delta variable
	    b(cons1) = 0;
	    %Second abs constaint:
	    N(cons1+1,nrxns+1) = d; %This is the normalization variable
	    N(cons1+1,k) = -1;  
            N(cons1+1,s2+1) = -1; %delta variable
	    b(cons1+1) = 0;

	    L(s2+1) = 0;         % this can be left as 0 in the CC transform 
	    U(s2+1) = infval;    % because it is just the same has having -delta <= 0
	    csense(cons1)   = 'L';
	    csense(cons1+1) = 'L';
	    f(s2+1) = - 1/s;
	    if m.rev(k)
	      N(cons1,k+1) = 1;
	      N(cons1+1,k+1) = -1;
	      k = k+1;
	    end
	    if k == r_group(k)
	      s1=s1+2;
	    end
	    s2=s2+1;
        end %end of if not nan
    end %end while k < nrxns
    disp('s1 s2:');
    disp([s1 s2]);
    %disp uptake value:
    %disp('Uptake 1:');
    %disp(m.lb(uptake(1)));

    %disp('s bounds');
    %disp([L(length(m.rxns)+1) U(length(m.rxns)+1)]);

    %disp('z bounds');
    %disp([L(length(m.rxns)+2) U(length(m.rxns)+2)]);

    disp(['Not Reversible: ' num2str(sum(~m.rev))]);
 %   if splitTurn
 %     tic; [v,fOpt,conv,vbasS,cbasS] = easyLP(f,N,b,L,U,csense,vbasS,cbasS); toc
 %   else
      tic; [v,fOpt,conv,vbasN,cbasN] = easyLP(f,N,b,L,U,csense,vbasN,cbasN); toc
 %   end
    %BB: We don't want to constrain the model to be at a particular flux,
    %but above we do initially normalize the expression to be on the same
    %scale as the flux, as well as setting an initial flux scale, so as
    %to help with the numerical computations taking place in the solver.
    %if firstRun
    %  for i = 1:length(gene_to_scale)
    %    uptake              = find(strcmp(gene_to_scale{i},m.rxnNames));
    %    m.lb(uptake) = upt_lb_save(i);
    %    m.ub(uptake) = upt_ub_save(i);
    %  end
    %  firstRun = false;
    %end
    corrval = fOpt;
    disp(fOpt);
    disp(v(nrxns+1));
    disp(v(nrxns+2));
    if conv
 %       if ~splitTurn
          v_orig = v;
          if v(nrxns+2) ~= 0
            v_orig = v/v(nrxns+2); %Transform to original
	  end
          v_sol = v_orig(1:nrxns);
	  v_all = [v_all v_sol];
   	  nvar = v_orig(nrxns+1);
	  %nvarp = nvar*nvarp;
	  %disp('New nvar, cumulative nvar is:');
	  %disp([nvar nvarp]);
	  disp('New nvar, zvar is:');
	  disp([nvar v(nrxns+2)]);
	  %Renormalize expression:
	  %r = r*nvar;
	  %r_sd = r_sd*nvar;
	%else
	  [m.lb m.ub m.rev] = setRxnDirection(v(1:nrxns), m.lb, m.ub, m.rev, nrxns, cnt);
	%end
      	disp(v_sol(1:15)');
        disp([sum(~m.rev) nR_old]);    
    end
%    splitTurn = ~splitTurn;
end % end of while sum(~m.rev) > nR_old


function [v,fOpt,conv,svbas,scbas] = easyLP(f,a,b,vlb,vub,csense,vbas,cbas)
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

% matlab can crash if inputs nan
if any(isnan(f))||any(any(isnan(a)))||any(isnan(b))...
        ||any(isnan(vlb))||any(isnan(vub))||any(isnan(csense)) 
    error('nan inputs not allowed');
end

% initialize
v = zeros(size(vlb));
v = v(:);
f = full(f(:));
vlb = vlb(:);
vub = vub(:);
% remove any tight contstraints as some solvers require volume > 0
j1 = (vlb ~= vub); 
j2 = (vlb == vub);

v(j2) = vlb(j2);
b = b(:) - a*v;
%b(isnan(b)) = 0; % likely a subtraction of nans
a(:,j2) = []; 
vlb(j2) = []; 
vub(j2) = [];
f0 = f;
f(j2) = [];
fOpt = nan;

if any(isnan(b))
    %disp(b);
    error('nan inputs not allowed: something went wrong');
end

params.method = 1;
if nargin > 6 && length(vbas) > 0
  %params.vbasis = zeros(1,length(vlb));
  %params.vbasis(1:length(vbas)) = vbas;
  params.vbasis = vbas;
  params.cbasis = zeros(1,length(b));
  params.cbasis(1:length(cbas)) = cbas;
end

solution = solveCobraLP(...
    struct('A',a,'b',b,'c',f,'lb',vlb,'ub',vub,'osense',-1,'csense',csense) , ...
     'GurobiParams',params);
%    'printLevel',1);

% define outputs
conv = solution.stat == 1;
svbas = []; %svbas = solution.basis;
scbas = []; %scbas = solution.cbasis;

if conv
    v0 = solution.full;
    v(j1) = v0;
    %v = v0;
    fOpt = f0'*v;
    %fOpt = f0'*v0;
    disp(['Convergent optimum is: ' num2str(solution.obj)]);
    if isnan(fOpt)
      disp('Converged, but fOpt still nan!');
    end
end



function [iLB iUB irev] = setRxnDirection(vI, iLB, iUB, irev, nrxns, cnt)
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

%An extremely small threshold doesn't help much: if the directionality is so uncertain based on 
%gene expression data, it should likely not be added to the objective function:

%rthresh = 0.5*exp(1/(cnt^2+0.5)); 
rthresh = 0.5;

%Another possibility is to incrementally lower this threshold

%nr_rxns = length(rev2irrev);
tol = 1e-9;


k = 0;
while k < nrxns 
  k = k+1;
  if irev(k)
    vSum = vI(k)+vI(k+1);
    if vI(k) >= tol && (vI(k) - vI(k+1) <= tol)
%    if 1 == 0
      iLB(k+1) = 0;
      iUB(k+1) = 0;
      irev(k) = 0;
      irev(k+1) = 0;

      iLB(k) = 0;
      iUB(k) = 0;
      irev(k) = 0;
      irev(k+1) = 0;      
    elseif vI(k) / vSum > rthresh
      %Forward reaction
      %rLB(i) = 0;
      iLB(k+1) = 0;
      iUB(k+1) = 0;

      %rrev(i) = 0;
      irev(k) = 0;
      irev(k+1) = 0;
    elseif vI(k+1) / vSum > rthresh
      %Backward reaction
      %rUB(i) = 0;
      iLB(k) = 0;
      iUB(k) = 0;

      %rrev(i) = 0;
      irev(k) = 0;
      irev(k+1) = 0;      
    end
    k = k+1;
  end
end