function newFilename = convertExpressionFilename(originalFilename)
newFilename=strrep(originalFilename,'(','_');
newFilename=strrep(newFilename,')','_');
newFilename=strrep(newFilename,' ','_');
newFilename=strrep(newFilename,'/','_');
newFilename=strrep(newFilename,'-','_');
end

