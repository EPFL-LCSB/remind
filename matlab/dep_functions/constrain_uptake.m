function model = constrain_uptake(model,rxns,rho,M)

%Adds the constrains on the uptake

%Get indexes

[~,idx_BU] = ismember(strcat('BU_',rxns),model.varNames);
[~,idx_MFA] = ismember(strcat('MFA_',rxns),model.varNames);
%idx_R = find(ismember(model.varNames,strcat('R_',rxns,'_reverse')));
[~,idx_R] = ismember(strcat('R_',rxns),model.varNames);
[~,idx_BFUSE] = ismember(strcat('BFUSE_',rxns),model.varNames);
[~,idx_FU] = ismember(strcat('FU_',rxns),model.varNames);

%idx_BU = getAllVar(model,{'BU'});

%% Coupling constraint MFA-BU<0

[NumCons, NumVars] = size(model.A);

for i=1:length(rxns)
    NewCons = zeros(NumVars,1);
    NewCons(idx_BU(i)) = -1;
    NewCons(idx_MFA(i)) = 1;
    model.A(NumCons+i,:) = NewCons;
    model.rhs(NumCons+i) = 0;
    model.constraintNames{NumCons+i} = strcat('Coupling_', num2str(i));
    model.constraintType{NumCons+i} = '<';
end

%% R + MFA(M-rho)<M

[NumCons, NumVars] = size(model.A);

for i=1:length(rxns)
    NewCons = zeros(NumVars,1);
    NewCons(idx_R(i)) = 1;
    NewCons(idx_MFA(i)) = M-rho;
    model.A(NumCons+i,:) = NewCons;
    model.rhs(NumCons+i) = M;
    model.constraintNames{NumCons+i} = strcat('MFA_Cut_', num2str(i));
    model.constraintType{NumCons+i} = '<';
end

%% BU -MFA +BFUSE +FU <1

[NumCons, NumVars] = size(model.A);

for i=1:length(rxns)
    NewCons = zeros(NumVars,1);
    NewCons(idx_BU(i)) = 1;
    NewCons(idx_MFA(i)) = -1;
    NewCons(idx_BFUSE(i)) = 1;
    NewCons(idx_FU(i)) = 1;
    model.A(NumCons+i,:) = NewCons;
    model.rhs(NumCons+i) = 1.5;
    model.constraintNames{NumCons+i} = strcat('Last_Cut_', num2str(i));
    model.constraintType{NumCons+i} = '<';
end
