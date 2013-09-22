function trstr = rstrtrim(instr, delim)
% Trims the character delim from the right end
% of a string.
    trstr = instr;
    while 1
        lstr = length(trstr);
	if lstr > 0
	    if trstr(lstr) == delim
		trstr = trstr(1:(lstr - 1));
                lstr = lstr - 1;
            else
                break;
	    end
	end
    end
   