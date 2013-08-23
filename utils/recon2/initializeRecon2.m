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

% Unconstrain reactions (namely boundary reactions) % This may be
% undesirable in some cases, but in such cases, more refined constraints
% should probably be used in any case.
%rec2.lb(modelIn.lb<0 & modelIn.lb>-1000) = -1000;
% Apparently not necessary in recon 2, but just in case:
%rec2.ub(modelIn.ub>0 & modelIn.ub<1000) = 1000;

%Need to set glucose uptake to ~ 3 mmol/gDWh