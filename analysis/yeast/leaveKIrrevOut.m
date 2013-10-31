function leaveKIrrevOut(model, modelOld, maxK, experiment)

corrDiffThresh = 0.03;
[modelAll, boundChanges] = useYN5irrevs(modelOld, model);

%We know that only lb's change.

N = length(boundChanges.LBOldIdx);
maxK = min(maxK, length(boundChanges.LBOldIdx));
corrOutFi = fopen(['leaveKOut_maxK_' num2str(maxK)], 'w');

%We need to record which rxns are removed to check for trend
%in correlations.

%But how do we analyze it? This is sort of like looking for 
%epistasis, where correlation is the fitness. 

%Another way is to assign to each rxn the pair (k, corr) when
%corr drops below a certain level, so that it can be easily plotted.
%But this ignores epistasis: (a,b,c) might given some plotted value for
% c even if (a,b,c) ~= (a,b). So we would need to further check this.

%But at this point, a recursive epistasis function is starting to make
% more sense. If the correlation drop is greater than 1%, record it, along
%with comma-separated list of old rxn indices and the new rxn index:
% r1,r2,r3:r4 \t corrVal - prevCorrVal \t corrVal \t prevCorrVal

%The downside to this is that it is more computation: all permutations
%instead of n choose k. Instead do n choose k

% r1,r2,r3 \t corrVal - prevCorrVal \t corrVal \t prevCorrVal

%To plot, x-axis: k, y-axis: corr diff, label: ??
%The other plotting method above would still work, but the insight
%here is that the label may be multiple reactions.

%In this case, the corr diff should probably be the maximal difference.
%Presumably, the non-maximal versions will still be plotted, but at
%smaller k (more significant) and smaller corr diffs (less signficant).

%Given this, we should use a Map Container with key num2str(kset(i,:))
%and value equal to the correlation.

corrMap = containers.Map();
for k = 1:maxK
    nCk = nchoosek(N, k);
    ksets = zeros(nCk, k);
    ksets(1, :) = [1:k];
    ksetsCorr = zeros(1, nCk);
    parfor i = 2:nCk
        ksets(i, :) = ksets(i - 1, :);
        ksets(i, k) = ksets(i, k) + 1; % necessary?
        cc = k;
        % First backtrack to the appropriate index and value
        while cc > 1
            if ksets(i, cc) > N + cc - k
                cc = cc - 1;
                ksets(i, cc) = ksets(i, cc) + 1;
                continue;
            else
                break;
            end
        end
        % Now continue counting up in the vector
        while cc < k
            cc = cc + 1;
            ksets(i, cc) = ksets(i, cc - 1) + 1;
        end

        %
        %Run simulation and get correlation.
        %
        [modelKout, ~] = useYN5irrevs(modelOld, model, ksets(i, :));
        corrCurr = yeastResults(modelKout, {'FALCON'}, experiment, false);
        ksetsCorr(i) = corrCurr;
    end
    for i = 2:nCk
        corrCurr = ksetsCorr(i);
        ksetiStr = num2str(ksets(i, :));
        corrMap(ksetiStr) = corrCurr;
        % Compare current correlation to 'ancestral' correlations
        maxCorrDiff = -inf;
        priorCorr = -inf;
        for rIdx = 1:k
            priorVec = ksets(i, :);
            priorVec(rIdx) = [];
            pvecStr = num2str(priorVec);
            if corrMap(pvecStr) - corrCurr > maxCorrDiff
                maxCorrDiff = corrMap(pvecStr) - corrCurr;
                priorCorr = corrMap(pvecStr);
            end
        end
        if maxCorrDiff > corrDiffThresh
            fprintf(corrOutFi, '%s\t%g\t%g\t%g\n', ksetiStr, ... 
                    maxCorrDiff, corrCurr, priorCorr);
        end
        i = i + 1;

    end % while i <= nCk
    
    %Debugging purposes:
    %dlmwrite(['ksets' num2str(k) '.csv'], boundChanges.LBOldIdx(ksets));
end
