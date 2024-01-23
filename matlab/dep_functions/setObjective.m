function [model_ILP] = setObjective(model,obj,type)
% function to set an objective function for the ILP part of ReMIND
% here competitiona and cooperation is implemented
model.f(:) = 0;

if strcmp(obj,'Cooperation')
    idx_xp = find(contains(model.varNames,'xp_'));
    model.f(idx_xp) = 1;
elseif strcmp(obj,'Competition')
    idx_xn = find(contains(model.varNames,'xn_'));
    model.f(idx_xn) = 1;
end

if strcmp(type,'max')
    model.objtype = -1;
elseif strcmp(type,'min')
    model.objtype = 1;
end

model_ILP = model;

end