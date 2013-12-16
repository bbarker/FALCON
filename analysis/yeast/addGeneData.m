function [n,n_sd] = addGeneData(g)

% kieran: 22 july 11

n = nan;
n_sd = nan;

ApmB = '[0-9\.]+±[0-9\.]+';

if ~isempty(g)
    while isnan(n)

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
