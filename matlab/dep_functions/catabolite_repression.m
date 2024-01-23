function model = catabolite_repression(model,CarbonSources,N)

%Adds constraints to BU of MainCarbonSources so that SUM BU<1 (so only one source is used) 

%Get idx of BU variables
[~,idx] = ismember(strcat('BU_',CarbonSources),model.varNames);
%Get size of constraint matrix A
[NumCons, NumVars] = size(model.A);
%Define new constraint
NewCons = zeros(NumVars,1);
NewCons(idx) = 1;
%Add constraint to constraint matrix
model.A(NumCons+1,:) = NewCons;
model.rhs(NumCons+1) = N+0.5;
model.constraintNames{NumCons+1} = 'catabolite_repr';
model.constraintType{NumCons+1} = '<';
