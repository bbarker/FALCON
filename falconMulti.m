function [v_sol, corrval, nvar, v_all, fTime, fIter, ...
    v_sol_s, corrval_s, nvar_s, fTime_s, fIter_s]  =  ...
    falconMulti(m, r, r_sd, r_group, nReps, rc, minFit,  ...
    EXPCON, FDEBUG, LPmeth)

% This is like falcon, but returns a mean for v_sol, corrval,
% nvar, fTime, and fIter.
% v_all is not intended to be used in this script, 
% so it is returned as [].
%
%INPUT    same as falcon.m, and also: 
%
% nReps   the number of replicate falcon fits.
%
%OUTPUT   same as falcom.m, and also:
% 
% v_sol_s, corrval_s, nvar_s, fTime_s, fIter_s
%    (these are the standard deviations of their namesakes)
%
% fileNameOut (falcon_nReps_StartTime.mat): distribution data for
%

nRxns = length(m.rxns);
v_sol_Dist   = zeros(nReps, nRxns);
corrval_Dist = zeros(nReps, 1);
nvar_Dist    = zeros(nReps, 1);
fTime_Dist   = zeros(nReps, 1);
fIter_Dist   = zeros(nReps, 1);
v_all = [];

timeInit = num2str(now());

%%%% Have to reproduce optional argument code:
% ?? Could this be included as a macro - just run a script?
if ~exist('rc', 'var')
    rc = 0;
end
if ~exist('minFit', 'var')
    minFit = 0;
end
if ~exist('FDEBUG', 'var')
    FDEBUG = false;
end
if ~exist('EXPCON', 'var')
    EXPCON = true;
end
if ~exist('LPmeth', 'var')
    LPmeth = 1; % dual-simplex for gurobi
end
%%%%

parfor i = 1:nReps
    [v_sol, corrval, nvar, ~, fTime, fIter] = ...
        falcon(m, r, r_sd, r_group, rc, minFit, EXPCON, FDEBUG, LPmeth);
    v_sol_Dist(i, :)  = columnVector(v_sol)';
    corrval_Dist(i)   = corrval;
    fTime_Dist(i)     = fTime;
    fIter_Dist(i)     = fIter;
end

v_sol   = mean(v_sol_Dist)';
nvar    = mean(nvar_Dist);
corrval = mean(corrval_Dist);
fTime   = mean(fTime_Dist);
fIter   = mean(fIter_Dist);

v_sol_s   = std(v_sol_Dist)';
nvar_s    = std(nvar_Dist);
corrval_s = std(corrval_Dist);
fTime_s   = std(fTime_Dist);
fIter_s   = std(fIter_Dist);

fileNameOut = ['falcon_' num2str(nReps) '_' timeInit '.mat'];
save(fileNameOut, 'v_sol', 'v_sol_s', 'v_sol_Dist', ...
    'corrval', 'corrval_s', 'corrval_Dist',         ...
    'nvar', 'nvar_s', 'nvar_Dist',                  ...
    'fTime', 'fTime_s', 'fTime_Dist',               ...
    'fIter', 'fIter_s', 'fIter_Dist');
