function [n,n_sd] = addGeneData(g)
% kieran: 22 july 11

n = nan;
n_sd = nan;

%ApmB = '[0-9\.]+±[0-9\.]+';
ApmB = '[0-9\.]+e?-?[0-9]*±[0-9\.]+e?-?[0-9]*';
loopcnt = 0;
geneNameLen = 5; %approximate lower bound on gene name length
if ~isempty(g)
    gLen = length(g);
    while (isnan(n) && loopcnt < floor(gLen/geneNameLen)) %arbitrary cutoff
        loopcnt = loopcnt + 1;
        %Checking for bad stuff, e.g. iMM1415
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
            match_expr      = '([0-9\.]+e?-?[0-9]*)±([0-9\.]+e?-?[0-9]*)';
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
            oldExp = regexp(g, match_expr, 'tokens', 'ignorecase');
            %replace_expr    = '${AandB($1,$2)}';
            if length(oldExp) >= 1
                replace_expr = str(eval(['AandB(' char(39) oldExp{1}{1} char(39) ...
                                  ',' char(39) oldExp{1}{2} char(39) ')']));
                g = regexprep(g,match_expr,replace_expr,'once','ignorecase');
            end
           
            % replace or
            match_expr      = ['(',ApmB,')\s+or\s+(',ApmB,')'];
            oldExp = regexp(g, match_expr, 'tokens', 'ignorecase');
            %replace_expr    = '${AorB($1,$2)}'
            if length(oldExp) >= 1
                replace_expr = str(eval(['AorB(' char(39) oldExp{1}{1} char(39) ...
                                  ',' char(39) oldExp{1}{2} char(39) ')']));
                g = regexprep(g,match_expr,replace_expr,'once','ignorecase');
            end
	    %format without brackets
            %if strcmp(g, g0)
            %    match_expr      = [ApmB,'\s+or\s+',ApmB];
            %    replace_expr    = '${AorB($1,$2)}';
            %    g = regexprep(g,match_expr,replace_expr,'once');
            %end
        end
    end
end

function s = str(s)
if isnumeric(s)
    s = num2str(s);
end    



function str = AandB(str1,str2) %#ok<DEFNU>

ApmB = '([0-9\.]+e?-?[0-9]*)±([0-9\.]+e?-?[0-9]*)';
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
ApmB = '([0-9\.]+e?-?[0-9]*)±([0-9\.]+e?-?[0-9]*)';
match_expr      = ApmB;
m1              = eval(regexprep(str1,match_expr,'$1'));
s1              = eval(regexprep(str1,match_expr,'$2'));
m2              = eval(regexprep(str2,match_expr,'$1'));
s2              = eval(regexprep(str2,match_expr,'$2'));

m = m1 + m2;

s = sqrt(s1^2 + s2^2);

str = [num2str(m),'±',num2str(s)];