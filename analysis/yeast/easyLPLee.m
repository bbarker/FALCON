function [v,fOpt,conv] = easyLPLee(f,a,b,vlb,vub)
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
        ||any(isnan(vlb))||any(isnan(vub))  
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
a(:,j2) = []; 
vlb(j2) = []; 
vub(j2) = [];
f0 = f;
f(j2) = [];
fOpt = nan;

% solve
csense(1:length(b)) = 'E';

solution = solveCobraLP(...
    struct('A',a,'b',b,'c',f,'lb',vlb,'ub',vub,'osense',-1,'csense',csense));

% define outputs
conv = solution.stat == 1;

if conv
    v0 = solution.full;
    v(j1) = v0;
    fOpt = f0'*v;
end
