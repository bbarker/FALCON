function words = strsplit(string, delimiter)
delimiter=sprintf(delimiter);
inds=strfind(string, delimiter);
words={};
if(length(inds)==0)
  words={string};
  return
end
for i=1:length(inds)
    if(i==1)
        words{end+1}=string(1:inds(i)-1);
    %elseif(i==length(inds))
        %words{end+1}=string(inds(i-1)+length(delimiter):end);
    else
        words{end+1}=string(inds(i-1)+length(delimiter):inds(i)-1);
    end
end
words{end+1}=string(inds(end)+1:end);
end

