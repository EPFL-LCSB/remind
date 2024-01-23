function [model_ILP] = changeConstraintParameter(model,P,new_value)
%the second input is the parameter to change, the third is the value to
%assign to that parameter.

if P == 'N'
    %change N
    model.N = new_value;
    %change yUse constraint
    idx_r = find(contains(model.constraintNames,'yUse_'));
    idx_c = find(contains(model.varNames,'y_'));
    for i=1:length(idx_r)
        [~,idx] = find(model.A(idx_r(i),idx_c));
        model.A(idx_r(i),idx_c(idx))=-new_value;
    end

elseif P == 'M'
    %change M
    model.M = new_value;

    %changeUpcomp1
    idx_r = find(contains(model.constraintNames,'UpComp1_'));
    idx_c = find(contains(model.varNames,'xn_'));
    for i=1:length(idx_r)
        [~,idx] = find(model.A(idx_r(i),idx_c));
        model.A(idx_r(i),idx_c(idx))=-new_value;
    end

    %changeUpcomp2
    idx_r = find(contains(model.constraintNames,'UpComp2_'));
    model.rhs(idx_r) = new_value - 2;
    for i=1:length(idx_r)
        [~,idx] = find(model.A(idx_r(i),idx_c));
        model.A(idx_r(i),idx_c(idx))=new_value;
    end

    %changeCoop2
    idx_r = find(contains(model.constraintNames,'Coop2_'));
    idx_c = find(contains(model.varNames,'mu_'));
    model.rhs(idx_r) = new_value;
    for i=1:length(idx_r)
        [~,idx] = find(model.A(idx_r(i),idx_c));
        model.A(idx_r(i),idx_c(idx))=new_value;
    end

    %changeCoop4
    idx_r = find(contains(model.constraintNames,'Coop4_'));
    idx_c = find(contains(model.varNames,'nu_'));
    model.rhs(idx_r) = new_value;
    for i=1:length(idx_r)
        [~,idx] = find(model.A(idx_r(i),idx_c));
        model.A(idx_r(i),idx_c(idx))=new_value;
    end

end

model_ILP = model;

end