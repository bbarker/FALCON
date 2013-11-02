function leaveKIrrevOut(model, modelOld, maxK, experiment)

corrDiffThresh = 0.05;
[modelAll, boundChanges] = useYN5irrevs(modelOld, model);

%We know that only lb's change.

N = length(boundChanges.LBOldIdx);
maxK = min(maxK, length(boundChanges.LBOldIdx));
corrOutFi = fopen([num2str(experiment) 'leaveKOut_maxK_' num2str(maxK)], 'w');

%We use the following output format for each line, where ri are reactions:
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
% Do 'WT' simulation (all K irrevs included):
corrMap('') = yeastResults(modelAll, {'FALCON'}, experiment, false);
WTfit = corrMap('');
for k = 1:maxK
    % The idea is that we start off with an initial ordered
    % vector and continue to the next ordered vector as we
    % go down the matrix, e.g., 1 2 3, 1 2 4, ... . This
    % creates all ordered sets of a given size k.
    nCk = nchoosek(N, k);
    ksets = zeros(nCk, k);
    ksets(1, :) = [1:k];
    ksetsCorr = zeros(1, nCk);
    for i = 2:nCk
        tmpVec = ksets(i - 1, :);
        %ksets(i, :) = ksets(i - 1, :);
        tmpVec(k) = tmpVec(k) + 1; % necessary?
        cc = k;
        % First backtrack to the appropriate index and value
        while cc > 1
            if tmpVec(cc) > N + cc - k
                cc = cc - 1;
                tmpVec(cc) = tmpVec(cc) + 1;
                continue;
            else
                break;
            end
        end
        % Now continue counting up in the vector
        while cc < k
            cc = cc + 1;
            tmpVec(cc) = tmpVec(cc - 1) + 1;
        end
        ksets(i, :) = tmpVec;
    end
    parfor i = 1:nCk
        %
        %Run simulation and get correlation.
        %
        leaveOut = boundChanges.LBOldIdx(ksets(i, :));
        [modelKout, ~] = useYN5irrevs(modelOld, model, leaveOut);
        corrCurr = yeastResults(modelKout, {'FALCON'}, experiment, false);
        ksetsCorr(i) = corrCurr;
        %ksetsCorr(i) = rand() %filler        
    end

    %Debugging purposes:
    dlmwrite(['ksets' num2str(k) '.csv'], ksets);
    %dlmwrite(['ksets' num2str(k) '.csv'], boundChanges.LBOldIdx(ksets));

    for i = 1:nCk
        corrCurr = ksetsCorr(i);
        ksetiStr = num2str(ksets(i, :));
        corrMap(ksetiStr) = corrCurr;
        % Compare current correlation to 'ancestral' correlations
        mut1Corr = -inf;
        mut2Corr = -inf;
        % not interested in 0 = self or 2^k-1 = WT
        if k > 1
            for psetIdx = 1:(2^k - 2)
                mut1RxnIdxs = ksets(i, :);
                mut2RxnIdxs = ksets(i, :);
                mut1SetIdxs = dec2bin(psetIdx, k);
                % One's complement:
                mut2SetIdxs = dec2bin(2^k - 1 - psetIdx, k);
                mut1RxnIdxs = mut1RxnIdxs(boolean(str2numvec(mut1SetIdxs)));
                mut2RxnIdxs = mut2RxnIdxs(boolean(str2numvec(mut2SetIdxs)));
                mut1Str = num2str(mut1RxnIdxs);
                mut2Str = num2str(mut2RxnIdxs);
                mut1Corr = corrMap(mut1Str);
                mut2Corr = corrMap(mut2Str);
                % If the negative epistasis is greater than a certain level:
                nEpistasis = mut1Corr * mut2Corr / (WTfit^2) ...
                             - corrCurr/WTfit;
                mut1RxnStr = num2str(boundChanges.LBOldIdx(mut1RxnIdxs));
                mut2RxnStr = num2str(boundChanges.LBOldIdx(mut2RxnIdxs));
                rxnIdxsStr = num2str(boundChanges.LBOldIdx(ksets(i, :)));
                if nEpistasis > corrDiffThresh
                    fprintf(corrOutFi, '%s\t%g\t%g\t%g\t%s\t%g\t%s\n', ...
                        rxnIdxsStr, nEpistasis, corrCurr, ...
                        mut1Corr, mut1RxnStr, mut2Corr, mut2RxnStr);
                end
            end
        else % separate case for single mutants
            for rIdx = 1:k
                corrCurr = ksetsCorr(i);
                ksetiStr = num2str(ksets(i, :));
                corrMap(ksetiStr) = corrCurr;
                fitDiff = 1 - corrCurr/WTfit;
                rxnIdxsStr = num2str(boundChanges.LBOldIdx(ksets(i, :)));
                if fitDiff > corrDiffThresh
                    fprintf(corrOutFi, '%s\t%g\t%g\t%g\t%s\t%g\t%s\n', ...
                        rxnIdxsStr, fitDiff, corrCurr, ...
                        WTfit, '', WTfit, '');
                end
            end
        end
        i = i + 1;
    end % while i <= nCk
end
