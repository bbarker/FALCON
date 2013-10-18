function [solutionDel, solutionWT] = linearMOMAOneShot(model, WTflux, objMult, ...
                                         fw, gw, rvec, rxnList, allowLoops)
%linearMOMA Performs a linear version of the MOMA (minimization of
%metabolic adjustment) approach 
%
% [solutionDel,solutionWT,totalFluxDiff,solStatus] = 
%       linearMOMA(modelWT,modelDel,osenseStr,minFluxFlag,verbFlab)
%
%INPUTS
% fw                regularization weight
% gw                LMOMA weight
% rvec              Perturbation weights for FBA objectives: 
%                   model.c(1:nRxns) += rvec
%
% modelWT           Wild type model
% modelDel          Deletion strain model
%
%OPTIONAL INPUTS
% osenseStr         Maximize ('max')/minimize ('min') (Default = 'max')
% minFluxFlag       Minimize the absolute value of fluxes in the optimal MOMA
%                   solution (Default = false)
% verbFlag          Verbose output (Default = false)
% 
%OUTPUTS
% solutionDel       Deletion solution structure
% solutionWT        Wild-type solution structure
% totalFluxDiff     Value of the linear MOMA objective, i.e. sum|v_wt-v_del|
% solStatus         Solution status
%
% Solves the problem
%
% min sum|v_wt - v_del|
%     S_wt*v_wt = 0
%     lb_wt <= v_wt <= ub_wt
%     c_wt'*v_wt = f_wt
%     S_del*v_del = 0
%     lb_del <= v_del <= ub_del
%
% Here f_wt is the optimal wild type objective value found by FBA
%
% Notes:
%
% 1) This formulation allows for selecting the most appropriate
% optimal wild type FBA solution as the starting point as opposed to
% picking an arbitrary starting point (original MOMA implementation).
%
% 2) The reaction sets in the two models do not have to be equal as long as
% there is at least one reaction in common
%
% Markus Herrgard 11/7/06

    % Minimize the absolute value of fluxes to 'avoid' loopy solutions
    % Solve secondary LP to minimize one-norm of |v|
    % Set up the optimization problem
    % f is the regularization parameter, g is the residual parameter
    % min -w*v_b + f*sum(delta) + g*sum(resid)
    % 1: S*v1 = 0
    % 3: delta+ >= -v1
    % 3: (new) v1 <= delta --> v1 - delta <= 0
    % 4: delta- >= v1
    % 4: (new) v1 >= -delta --> v1 + delta >= 0
    % x5x (don't need for MOMA): c'v1 >= f (optimal value of objective)
    % 5: resid <= t1 --> v1 - b <= t1 --> v1 - t1 <= b
    % 6: resid >= -t1 --> v1 - b >= -t1 --> v1 + t1 >= b
    % delta+,delta- >= 0
    %One potential issue with this formulation is that it adds more contstraints
    %though it adds variables (only one "delta" instead of 2, etc.).

nMets = length(model.mets);
nRxns = length(model.rxns);

%Remember to remove/change this when needed
%fw = log(log(objMult)+1)/100;
%
if nargin < 8
  allowLoops = true;
end
    if (~isfield(model,'b'))
      LPproblem.b = zeros(size(model.S,1),1);
    else
      LPproblem.b = model.b;
    end

    LPproblem2.A = [model.S sparse(nMets,2*nRxns);
      speye(nRxns,nRxns) -speye(nRxns,nRxns) sparse(nRxns,nRxns);
      speye(nRxns,nRxns) speye(nRxns,nRxns) sparse(nRxns,nRxns);
      speye(nRxns,nRxns) sparse(nRxns,nRxns) -speye(nRxns,nRxns);
      speye(nRxns,nRxns) sparse(nRxns,nRxns) speye(nRxns,nRxns) ];
    LPproblem2.c = [zeros(nRxns,1);fw*ones(nRxns,1);gw*ones(nRxns,1)];
    %do not weight reactions that aren't part of rxnList:
    %not efficient, but simple:
    if nargin > 6
      LPproblem2.c(2*nRxns+1:end) = zeros(nRxns,1);
      for i = 1:length(rxnList)
        LPproblem2.c(2*nRxns + rxnList(i)) = gw;
      end
    end
    %biomass/growth correction:
    LPproblem2.c(nRxns+find(model.c))=0;
    LPproblem2.c(2*nRxns+find(model.c))=0;
    if nargin > 5
      LPproblem2.c(1:nRxns) = LPproblem2.c(1:nRxns) + columnVector(rvec);
    end
    LPproblem2.c(find(model.c)) = -objMult;
    LPproblem2.lb = [model.lb;zeros(2*nRxns,1)];
    LPproblem2.ub = [model.ub;10000*ones(2*nRxns,1)];
    LPproblem2.b = [LPproblem.b;zeros(2*nRxns,1);columnVector(WTflux);columnVector(WTflux)];
    if ~isfield(model,'csense')
        % If csense is not declared in the model, assume that all
        % constraints are equalities.
        LPproblem2.csense(1:nMets) = 'E';
    else % if csense is in the model, move it to the lp problem structure
        if length(model.csense)~=nMets,
            warning('Length of csense is invalid! Defaulting to equality constraints.')
            LPproblem2.csense(1:nMets) = 'E';
        else
            LPproblem2.csense = columnVector(model.csense);
        end
    end
    LPproblem2.csense((nMets+1):(nMets+nRxns)) = 'L';
    LPproblem2.csense((nMets+nRxns+1):(nMets+2*nRxns)) = 'G'; 
    LPproblem2.csense((nMets+2*nRxns+1):(nMets+3*nRxns)) = 'L'; 
    LPproblem2.csense((nMets+3*nRxns+1):(nMets+4*nRxns)) = 'G'; 
    LPproblem2.csense = columnVector(LPproblem2.csense);
    LPproblem2.osense = 1;
    if allowLoops
        solution = solveCobraLP(LPproblem2); %,printLevel,minNorm);
        solution.dual = []; % slacks and duals will not be valid for this computation. %BB:?
        solution.rcost = [];
    else
        MILPproblem2 = addLoopLawConstraints(LPproblem, model, 1:nRxns);
        solution = solveCobraMILP(MILPproblem2);
    end
    if (solution.stat > 0)
      solutionDel.x = solution.full(1:nRxns);
      solutionDel.f = sum(model.c.*solutionDel.x); 
    else
      solutionDel.x = WTflux*0;
      solutionDel.f = 0;  
    end

solutionWT.x = WTflux;
solutionWT.f = model.c'*WTflux; %Assumes same objective for both models.
solutionWT.stat = 1;
solutionWT.stat = solution.stat;
solutionDel.stat = solution.stat;
solStatus = solution.stat;
