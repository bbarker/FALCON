function [inactiveRxns,time,numOpts] = checkModelConsistency(model,r,C,fast)

% This function is designed to quickly identify dead-end reactions in a
% stoichiometric model. The algorithm is largely based on the heuristic
% speed-up to Flux Variability Analysis (FVA) proposed by Jerby et al. [1],
% with modifications to further reduce computation time in Matlab. The
% function can operate independently to report the inactive reactions for
% an entire model, or within a pruning algorithm (e.g., MBA) to examine the
% effect of removing reactions.

% INPUT
%  model                COBRA model structure
%  r                    name of reaction to be removed (for model pruning)
%  C                    list of core reaction names (for model pruning)

% When using the function independently, no additional inputs need to be
% specified besides the model.

if nargin < 2
    r = [];
    C = {};
end
if nargin < 4
    fast = 0;
end
if numel(r)
   % Remove reaction r from the model
    model = removeRxns(model,r); 
end
model.c(logical(model.c)) = 0;

inactiveRxns = r;

numOpts = 0;
t0 = clock;

% First check whether any core reactions are blocked by the removal of r.
rxnList = C;

% Checking for metabolite dead ends is accomplished entirely by matrix
% operations, and is therefore very fast in Matlab. If any core reaction
% contains a dead-end metabolite, the reaction itself will be a dead end.
% This check potentially avoids sequential optimizations, as the function
% can exit if any blocked core reactions are detected.
deadEndMets = detectDeadEnds_fast(model);
deadEnd = sum(full(model.S(deadEndMets,:)~=0),1)>0;
deadEndRxns = model.rxns(deadEnd);
deadEnd_C = intersect(rxnList,deadEndRxns);

% If no core reactions were found to be blocked based on metabolite dead
% ends, maximize and minimize reactions to identify those with zero flux
% capacity. If the option is specified, fastFVA is used to quickly scan
% through all reactions. 
if numel(deadEnd_C)
    inactiveRxns = union(deadEnd_C,inactiveRxns);
elseif fast
    display('Checking all reactions (fastFVA)...')
    model.c(logical(model.c)) = 0;
    [optMin,optMax] = fastFVA(model,0,'max','cplex');
    fastInactive = (abs(optMax) < 1e-6) & (abs(optMin) < 1e-6);
    fastInactiveRxns = model.rxns(fastInactive);
    inactiveRxns = union(inactiveRxns,fastInactiveRxns);
    
% Othwerise, this step follows the heuristic speed-up to FVA proposed by
% Jerby et al. to scan through core reactions only; if any core reactions
% are blocked, the function can exit.
else
    display('Checking core reactions...')
    changeCobraSolver('glpk');
    while numel(rxnList)
        numRxnList = numel(rxnList);
        model = changeObjective(model,rxnList);

        % Maximize all
        FBAsolution = optimizeCbModel(model,'max');
        numOpts = numOpts + 1;
        optMax = FBAsolution.x;
        active = (abs(optMax) >= 1e-6);
        activeRxns = model.rxns(active);
        rxnList = setdiff(rxnList,activeRxns);

        % Minimize all
        if ~numel(activeRxns)
            FBAsolution = optimizeCbModel(model,'min');
            numOpts = numOpts + 1;
            optMin = FBAsolution.x;
            active = (abs(optMin) >= 1e-6);
            activeRxns = model.rxns(active);
            rxnList = setdiff(rxnList,activeRxns);
        end

        numRemoved = numRxnList - numel(rxnList);

        if ~numRemoved
            randInd = randperm(numel(rxnList));
            i = rxnList(randInd(1));
            model = changeObjective(model,i);

            % Maximize reaction i
            FBAsolution = optimizeCbModel(model,'max');
            numOpts = numOpts + 1;
            optMax = FBAsolution.f;

            % Minimize reaction i (if reversible)
            lb_i = model.lb(strmatch(i,model.rxns,'exact'));
            if lb_i < 0 && optMax < 1e-6
                FBAsolution = optimizeCbModel(model,'min');
                numOpts = numOpts + 1;
                optMin = FBAsolution.f;
            else optMin = 0;
            end

            if (abs(optMax) < 1e-6) && (abs(optMin) < 1e-6)
                inactiveRxns = union(inactiveRxns,i);
                break;
            end

            rxnList = setdiff(rxnList,i);
        end
    end
end

% If all core actions remain active, r can be removed from the model. This
% step checks all other non-core reactions to determine whether anything
% else can be concurrently removed. Note: core reactions are included in
% the objective function, as this seems to speed up the algorithm, but are
% not individuall maximized or minimized (since this is redundant witht he
% above step).
if numel(inactiveRxns) <= 1 && fast == 0  
    display('Checking non-core reactions...')
    rxnList = model.rxns;
    inactiveRxns = union(inactiveRxns,deadEndRxns);
    rxnList = setdiff(rxnList,inactiveRxns);

    while numel(rxnList)
        numRxnList = numel(rxnList);
        model = changeObjective(model,rxnList);

        % Maximize all
        FBAsolution = optimizeCbModel(model,'max');
        numOpts = numOpts + 1;
        optMax = FBAsolution.x;
        active = (abs(optMax) >= 1e-6);
        activeRxns = model.rxns(active);
        rxnList = setdiff(rxnList,activeRxns);

        % Minimize all
        if ~numel(activeRxns)
            FBAsolution = optimizeCbModel(model,'min');
            numOpts = numOpts + 1;
            optMin = FBAsolution.x;
            active = (abs(optMin) >= 1e-6);
            activeRxns = model.rxns(active);
            rxnList = setdiff(rxnList,activeRxns);
        end

        numRemoved = numRxnList - numel(rxnList);

        if ~numRemoved
            randInd = randperm(numel(rxnList));
            i = rxnList(randInd(1));

            if ~ismember(i,C)
                model = changeObjective(model,i);

                % Maximize reaction i
                FBAsolution = optimizeCbModel(model,'max');
                numOpts = numOpts + 1;
                optMax = FBAsolution.f;

                % Minimize reaction i (if reversible)
                lb_i = model.lb(strmatch(i,model.rxns,'exact'));
                if lb_i < 0 && optMax < 1e-6
                    FBAsolution = optimizeCbModel(model,'min');
                    numOpts = numOpts + 1;
                    optMin = FBAsolution.f;
                else optMin = 0;
                end

                if (abs(optMax) < 1e-6) && (abs(optMin) < 1e-6)
                    inactiveRxns = union(inactiveRxns,i);
                end
            end

            rxnList = setdiff(rxnList,i);
            if ~numel(setdiff(rxnList,C))
                break;
            end
        end
    end
end

time = etime(clock,t0);
display(['checkModelConsistency time: ',num2str(time,'%1.2f'),' s'])