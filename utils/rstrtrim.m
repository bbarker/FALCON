function trstr = rstrtrim(instr, delim)
% Trims the character delim from the right end
% of a string.
    trstr = instr;
    lstr = length(trstr);
    while 1
	if lstr > 0
	    if instr(lstr) == delim
		trstr = trstr(1:(lstr - 1));
                lstr = lstr - 1;
            else
                break;
	    end
	end
    end
   