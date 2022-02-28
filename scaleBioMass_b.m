function model = scaleBioMass_b(model,component,new_value,balance_out)
  % scaleBioMass
  %   Scales the biomass composition
  %
  %   model          (struct) metabolic model in COBRA format
  %   component      (str) name of the component to rescale (e.g. "protein")
  %   new_value      (float) new total fraction for said component
  %   balance_out    (str, opt) if chosen, the name of another component with which
  %                  the model will be balanced out so that the total mass remains = 1 g/gDW
  %
  %   model          (struct) modified model
  %
  %   Usage: model = scaleBioMass(model,component,new_value,balance_out)
  %

%Measure current composition and rescale:
[~,P,C,R,D,L,I,F] = sumBioMass(model);
content_all = {'carbohydrate','protein','lipid','RNA','DNA','ion','cofactor'};
content_Cap = {'C','P','L','R','D','I','F'};
pos         = strcmp(content_all,component);
old_value   = eval(content_Cap{pos});
f           = new_value / old_value;
model       = rescalePseudoReaction(model,component,f);

%Balance out (if desired):
if nargin > 3
    diff_names = setdiff(setdiff(content_all,component), {'ion','cofactor'});
    diff = (new_value - old_value)/length(diff_names);
    for i = 1:length(diff_names)
        pos           = strcmp(content_all,diff_names{i});
        balance_value = eval(content_Cap{pos});
        f             = (balance_value - diff) / balance_value;
        model         = rescalePseudoReaction(model,diff_names{i},f);
    end
end

end
