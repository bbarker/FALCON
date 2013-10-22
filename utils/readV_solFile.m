function [v_solex v_solrev v_solirrev] = readV_solFile(fluxFileLoc)
    inputFI = fopen(fluxFileLoc);
    line=fgetl(inputFI);
    v_solex=[];v_solrev=[];v_solirrev=[];
    flagexfluxes=0;flagrevfluxes=0;flagirrevfluxes=0;
    while line~=-1
        if(~isempty(regexp(line,'v_sol')))
 	    if(~isempty(regexp(line,'v_solex')))
                flagexfluxes=1;flagrevfluxes=0;flagirrevfluxes=0;
	    elseif(~isempty(regexp(line,'v_solrev')))
	        flagexfluxes=0;flagrevfluxes=1;flagirrevfluxes=0;
	    elseif(~isempty(regexp(line,'v_solirrev')))
	      flagexfluxes=0;flagrevfluxes=0;flagirrevfluxes=1;
	    end
	else
	    if(flagexfluxes)
	      startindex=regexp(line,'(\-|\d|\.)+$');
            v_solex=[v_solex; str2num(line(startindex:length(line)))];
	    elseif(flagrevfluxes)
	      startindex=regexp(line,'(\-|\d|\.)+$');
            v_solrev=[v_solrev; str2num(line(startindex:length(line)))];
	    elseif(flagirrevfluxes)
	      startindex=regexp(line,'(\-|\d|\.)+$');
            v_solirrev=[v_solirrev; str2num(line(startindex:length(line)))];
	     end
        end
        line=fgetl(inputFI);
    end
    fclose(inputFI);
end