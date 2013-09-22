function newFilename = convertExpressionFileName(originalFilename)
newFilename = strrep(originalFilename, '(', '_');
newFilename = strrep(newFilename, ')', '_');
newFilename = strrep(newFilename, ' ', '_');
newFilename = strrep(newFilename, '/', '_');
newFilename = strrep(newFilename, '-', '_');
end

