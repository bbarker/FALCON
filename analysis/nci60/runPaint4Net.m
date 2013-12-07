inputFI=fopen('../rxnsToShow.txt');
line=fgetl(inputFI);
rxnsToShow={};
k=0;
while line~=-1
  k=k+1;
  rxnsToShow{k}=line;
  line=fgetl(inputFI);
end
draw_by_rxn(model,rxnsToShow,'true');