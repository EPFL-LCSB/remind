function model=block_essential(model,EssList)
%Adds a constraint to the sources in the essential list that doesn't allow BU 0 and MFA 1 -> so doesn't allow big uptake of essentials

%Get size of A
[NumCons, NumVars] = size(model.A);

% Add constraint -BU+MFA > 0

    for i=1:length(EssList)
        idx_BU = find(ismember(model.varNames,strcat('BU_',EssList(i))));
        idx_MFA = find(ismember(model.varNames,strcat('MFA_',EssList(i))));
        NewCons = zeros(NumVars,1);
        NewCons(idx_BU) = -1;
        NewCons(idx_MFA) = 1;
        model.A(NumCons+i,:) = NewCons;
        model.rhs(NumCons+i) = 0;
        model.constraintNames{NumCons+i} = strcat('Ess_Cut_', num2str(i));
        model.constraintType{NumCons+i} = '>';
    end