function [v_sol, lIter] = dataToFluxFix(m,r,r_sd)

% kieran: 21 sep 11

rev = false(size(m.rxns));

nR_old = 0;

% m.lb(~ismember(m.lb,[-inf,0,inf,-1000,1000])) = -inf;
% m.ub(~ismember(m.ub,[-inf,0,inf,-1000,1000])) = inf;

v_sol = zeros(size(m.rxns));

nRxns = length(m.rxns);

lIter = 0;
while sum(~m.rev) > nR_old
    lIter = lIter + 1;
    nR_old = sum(~m.rev);
    
    % 1. fit to data
    
    N = m.S;    
    L = m.lb;
    U = m.ub;
    f = zeros(size(m.rxns))';
    b = zeros(size(m.mets));
    
    kPerm = randperm(nRxns);
    for k = 1:nRxns
        kp = kPerm(k);
        %kp = k; %use for "deterministic" method.
        d = r(kp);
        s = r_sd(kp);
        if ~m.rev(kp) && ~isnan(d) && s>0
            [s1,s2] = size(N);
            N(s1+1, kp) = 1; N(s1+1, s2+1) = -1; N(s1+1, s2+2) = 1;
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
                    
                    %if conv && abs(fOpt)<=0 % irreversibly backward
                    if conv && (fOpt <= 0)    
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


