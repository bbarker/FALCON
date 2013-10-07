function testTestModels()
% Runs some simple tests on models as a sanity check
% for the models themselves. For descriptions of
% individual models see the appropriate model-generating
% .m file.


disp('*** Testing InTriangleOut model ***');

m = makeTestModel_InTriangleOut();
sol = optimizeCbModel(m, 'max', 'one');

F_3 = find(strcmp(m.rxns, 'F_3'));
F_4 = find(strcmp(m.rxns, 'F_4'));

E_I = find(strcmp(m.rxns, 'E_I'));
E_J = find(strcmp(m.rxns, 'E_J'));
E_C = find(strcmp(m.rxns, 'E_C'));

if sol.x(F_3) == 0 && sol.x(F_4) == 0
    disp('Test succeeded for no futile cycles in min 1-norm FBA');
else
    disp('Test FAILED for no futile cycles in min 1-norm FBA'); 
end

if sol.x(E_C) == (1/2) * sol.x(E_I)
    disp('Test succeeded for input-output stoichiometry in FBA');
else
    disp('Test FAILED for input-output stoichiometry in FBA');
end
disp(' ');
disp(' ');