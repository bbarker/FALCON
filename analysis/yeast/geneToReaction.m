function [r,r_sd,rxn_missing_gene] = geneToReaction(m,g,t,t_sd)

% kieran: 16 sep 11

% brandon 07 aug 12

r       = zeros(size(m.rxns));
r_sd    = zeros(size(m.rxns));
rxn_missing_gene = zeros(size(m.rxns));
true_missing = 0;
for k = 1:length(g)
    g{k} = strrep(g{k},'-','_');
end

for k = 1:length(m.rxns)
    ga = m.grRules{k};
    ga = strrep(ga,'-','_');    
    w = regexp(ga,'\<\w*\>','match'); 
    w = setdiff(w,{'and','or','AND','OR'});
    for kk = 1:length(w)
	j = find(strcmp(w{kk},g));
	if length(j) == 0 %for names encloses in parentheses:
	  j = find(strcmp(['(',w{kk},')'],g));
	end
	if numel(j) > 1
	  j = j(1); %temporary fix
	end
        n = t(j);
        n_sd = t_sd(j);
	if ~isempty(n)
          ga = regexprep(ga,['\<',w{kk},'\>'],[num2str(n),'±',num2str(n_sd)]); % ±
	  %Now try with parentheses ...
          ga = regexprep(ga,['\<','\(',w{kk},'\)','\>'],[num2str(n),'±',num2str(n_sd)]); % ±
	else
	  true_missing = true_missing+1;
	  %try right first
	  gatmp = regexprep(ga,['\<',w{kk},'\>','\s+and\s+'], '', 'ignorecase');
	  gatmp = regexprep(gatmp,['\<',w{kk},'\>','\s+or\s+'], '', 'ignorecase');

	  %try left
	  gatmp = regexprep(gatmp,['\s+and\s+','\<',w{kk},'\>'], '', 'ignorecase');
	  gatmp = regexprep(gatmp,['\s+or\s+','\<',w{kk},'\>'], '', 'ignorecase');

	  %Now try with parentheses ...
	  %try right first
	  gatmp = regexprep(gatmp,['\<','\(',w{kk},'\)','\>','\s+and\s+'], '', 'ignorecase');
	  gatmp = regexprep(gatmp,['\<','\(',w{kk},'\)','\>','\s+or\s+'], '', 'ignorecase');

	  %try left
	  gatmp = regexprep(gatmp,['\s+and\s+','\<','\(',w{kk},'\)','\>'], '', 'ignorecase');
	  gatmp = regexprep(gatmp,['\s+or\s+','\<','\(',w{kk},'\)','\>'], '', 'ignorecase');

	  if strcmp(gatmp,ga)
  	    rxn_missing_gene(k) = 1;
            disp(['skipping rxn# ' num2str(k) ' due to ' gatmp]);
	  else
	    ga = gatmp;
	  end
	end
	%disp(ga);
    end
    if ~(rxn_missing_gene(k))
      [n,n_sd] = addGeneData(ga);
      r(k) = n;
      r_sd(k) = n_sd;
    end
end
disp(['Number of reactions actually missing data: ', num2str(true_missing)]);
