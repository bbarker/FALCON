function trstr = lstrtrim(instr, delim)
% Trims the character delim from the left end
% of a string.
    trstr = instr;
    while 1
        lstr = length(trstr);
	if lstr > 0
	    if trstr(1) == delim
		trstr = trstr(2:lstr);
                lstr = lstr - 1;
            else
                break;
	    end
	end
    end
   