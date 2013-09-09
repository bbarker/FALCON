function [trimmedX trimmedY]=trimOutliers(X,Y)
[sortedX sortedXIdxs]=sort(X);
[sortedY sortedYIdxs]=sort(Y);
idxsToKeep=intersect(sortedXIdxs(3:end-3),sortedYIdxs(3:end-3));
trimmedX=X(idxsToKeep);
trimmedY=Y(idxsToKeep);
end