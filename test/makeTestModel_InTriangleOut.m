function model = makeTestModel_InTriangleOut()
%
% This is the model presented in Figure 1 of
% Computational approaches for understanding energy metabolism.
% See https://app.box.com/s/4xrz35so0s59v2z2l2r4
%
% This model's gene rules totally ignores complicated 
% gene rule interactions. Instead, the purpose is to test
% the optimization routine.

model.description = 'InTriangleOut';

model.S = [
1	-1	0	1	0	0	0	0	0	0;
0	0	1	-1	0	0	0	0	0	0;
0	1	-1	0	0	0	-1	0	0	0;
-2	0	0	0	1	0	0	0	0	0;
-1	0	0	0	0	1	0	0	0	0;
0	0	0	0	-1	0	0	1	0	0;
0	0	0	0	0	-1	0	0	1	0;
0	0	0	0	0	0	1	0	0	-1
];

[nmets, nrxns] = size(model.S);

model.rxns = {               ...
'F_1', 'F_2', 'F_3', 'F_4',  ...
'E_I', 'E_J', 'E_C',         ...
'B_I', 'B_J', 'B_C'          }';
model.rxnNames = model.rxns;

% Assume all reactions are reversible with one exception, 
% unlike in the figure illustration.
model.lb = -1*ones(nrxns, 1);
model.ub = ones(nrxns, 1);
model.rev = ones(nrxns, 1);
% The exception will define an uptake, 'E_I'
E_I = find(strcmp(model.rxns, 'E_I'));
model.rev(E_I) = 0;
model.lb(E_I) = 0;

% Default objective, flux out through the boundary 
% (B_c in the figure)
model.c = zeros(nrxns, 1);
model.c(find(strcmp(model.rxns, 'B_C'))) = 1;


model.mets = {                      ...
'A_c', 'B_c', 'C_c', 'I_c', 'J_c',  ...
'I_e', 'J_e', 'C_e'                 }';
model.metNames = model.mets;

% Now we add some default gene rules that aren't
% included in the paper - but only for internal
% reactions. 
%
% We use the following default scenario: 
% 1) All internal reactions have a single gene.
% 2) E_I is the only transporter requiring a gene.
model.grRules = {                ...
'GF_1', 'GF_2', 'GF_3', 'GF_4',  ...
'GE_I', '', '',                  ...
'', '', ''                       }';
model.rules = {                  ...
'x(1)', 'x(2)', 'x(3)', 'x(4)',  ...
'x(5)', '', '',                  ...
'', '', ''                       }';


% A simple method: won't handle complexes/isozymes:
model.genes = model.grRules(find(cellfun(@length, model.grRules)));
ngenes = length(model.genes);
model.rxnGeneMat = zeros(nrxns, ngenes);
for i = 1:ngenes
    model.rxnGeneMat(i, i) = 1;
end

