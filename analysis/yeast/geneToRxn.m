function [rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToRxn(model,genedata_filename);

% load transcript data
genedata	= importdata(genedata_filename);
genenames	= genedata.textdata(:,1);
genenames(1)= [];
gene_exp	= genedata.data(:,1);
gene_exp_sd	= genedata.data(:,2);

% map gene weighting to reaction weighting
[rxn_exp,rxn_exp_sd,rxn_missing_gene] = geneToReaction(model,genenames,gene_exp,gene_exp_sd);
% sds 0 -> small
rxn_exp_sd(rxn_exp_sd == 0) = min(rxn_exp_sd(rxn_exp_sd>0))/2;


%Some of the symbols used below (in geneToReaction
%and addGeneData) are not copied correctly in emacs.
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
	    disp(gatmp);
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
    else
      %disp('skipping rxn');
    end
end
disp(['Number of reactions actually missing data: ', num2str(true_missing)]);

function [n,n_sd] = addGeneData(g)
% kieran: 22 july 11

n = nan;
n_sd = nan;

ApmB = '[0-9\.]+±[0-9\.]+';

if ~isempty(g)
    while isnan(n)
        %Checking for bad stuff, cough cough iMM1415 cough
        gtmp = regexprep(g,'\s+or\s+or\s+', ' or ', 'ignorecase');
        gtmp = regexprep(gtmp,'\s+and\s+and\s+', ' and ', 'ignorecase');
	if strcmp(gtmp,g) == 0
	  disp('Warning: adjacent logical operators of the same type replaced with a single operator:');
	  disp(g);
	end
	g = gtmp;
        gtmp = regexprep(g,'and\s*$', '', 'ignorecase');
        gtmp = regexprep(gtmp,'^\s*and', '', 'ignorecase');
        gtmp = regexprep(gtmp,'or\s*$', '', 'ignorecase');
        gtmp = regexprep(gtmp,'^\s*or', '', 'ignorecase');
	if strcmp(gtmp,g) == 0
	  disp('Warning: trailing logical operators found and removed:');
	  disp(g);
	end
	g = gtmp;

        try 
            match_expr      = '([0-9\.])+±([0-9\.]+)';
            g_av            = regexprep(g,match_expr,'$1');
            g_sd            = regexprep(g,match_expr,'$2');
            n               = eval(g_av);
            n_sd            = eval(g_sd);
        catch %#ok<CTCH>
            % replace brackets
            match_expr      = ['\((\s*',ApmB,'\s*)\)'];
            replace_expr    = '$1';
            g = regexprep(g,match_expr,replace_expr);
            % replace and
            match_expr      = ['(',ApmB,')\s+and\s+(',ApmB,')'];
            replace_expr    = '${AandB($1,$2)}';
            g = regexprep(g,match_expr,replace_expr,'once','ignorecase');

            % replace or
            match_expr      = ['(',ApmB,')\s+or\s+(',ApmB,')'];
            replace_expr    = '${AorB($1,$2)}';
            g = regexprep(g,match_expr,replace_expr,'once','ignorecase');
	    %format without brackets
            %match_expr      = [ApmB,'\s+or\s+',ApmB];
            %replace_expr    = '${AorB($1,$2)}';
            %g = regexprep(g,match_expr,replace_expr,'once');
           
        end
    end
end

function str = AandB(str1,str2) %#ok<DEFNU>

ApmB = '([0-9\.])+±([0-9\.]+)';
match_expr      = ApmB;
m1              = eval(regexprep(str1,match_expr,'$1'));
s1              = eval(regexprep(str1,match_expr,'$2'));
m2              = eval(regexprep(str2,match_expr,'$1'));
s2              = eval(regexprep(str2,match_expr,'$2'));

[m,j] = min([m1,m2]);

if j == 1
    s = s1;
else
    s = s2;
end

str = [num2str(m),'±',num2str(s)];

function str = AorB(str1,str2) %#ok<DEFNU>
ApmB = '([0-9\.])+±([0-9\.]+)';
match_expr      = ApmB;
m1              = eval(regexprep(str1,match_expr,'$1'));
s1              = eval(regexprep(str1,match_expr,'$2'));
m2              = eval(regexprep(str2,match_expr,'$1'));
s2              = eval(regexprep(str2,match_expr,'$2'));

m = m1 + m2;

s = sqrt(s1^2 + s2^2);

str = [num2str(m),'±',num2str(s)];