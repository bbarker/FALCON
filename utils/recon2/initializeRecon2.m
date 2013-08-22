function rec2 = initializeRecon2(modelIn) 
% This script is used to set up default constraints
% given any version of Human Recon 2.XX.

rec2 = modelIn;

% Remove unhelpful transcript labels
for i = 1:length(rec2.rxns)
  rec2.grRules{i} = regexprep(modelIn.grRules{i},'(\d+)\.\d*','$1');
end
for i = 1:length(rec2.genes)
  rec2.genes{i} = regexprep(modelIn.genes{i},'(\d+)\.\d*','$1');
end

%Assign description if it doesn't exist
if ~isfield(rec2,'description')
  rec2.description = 'human_recon_2'
end

%Unconstrain reactions (namely boundary reactions)