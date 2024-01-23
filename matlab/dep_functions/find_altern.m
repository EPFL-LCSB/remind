function [model,main_c,DPs]=find_altern(model,MaxAlt,idx_BFUSE,idx_BU,idx_MFA)

sol = solveTFAmodelCplex(model);

main_c = {};

alt = 1;

Var_Names = model.varNames;

DPs = cell2table(Var_Names);

%min_size_c = length(intersect(idx_BU(sol.x(idx_BFUSE)<0.1),idx_BU(sol.x(idx_BU)>0.9)));

while ((alt <= MaxAlt) && ~(isempty(sol.x)) && ~(sol.val==0) )

    %Save solution

    name= strcat('Alternative_', num2str(alt));

    DPs.(name) = sol.x;

    fprintf('Number of alts:\t%d\n',alt);

    %Add integer cut

    [NumCons, NumVars] = size(model.A);

    actuse=intersect(idx_BU(sol.x(idx_MFA)<0.1),idx_BU(sol.x(idx_BU)>0.9));

    %if (abs(length(actuse)-min_size_c) > StopCriteria)
        %break
    %end

    NewCons = zeros(NumVars,1);
    NewCons(actuse) = 1;

    %add constraint to constraint matrix
    model.A(NumCons+1,:) = NewCons;
    model.rhs(NumCons+1) = 0.5;
    %length(actuse)-0.5;
    model.constraintNames{NumCons+1} = strcat('int_cut_', num2str(alt));
    model.constraintType{NumCons+1} = '<';
  
    main_c = [main_c; strrep(model.varNames(actuse),'BU_','')];

    main_c = unique(main_c,'stable');

    if length(actuse)<1
        break
    end

    sol = solveTFAmodelCplex(model);
        
        if isempty(sol.x) || (sol.val==0)
            break
        end

    alt=alt+1;
end

main_c = unique(main_c,'stable');
