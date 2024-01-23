function CarbonAtoms = getCarbonAtoms(model,metsIdx)

for i = 1:length(metsIdx)
    formula = char(model.metFormulas(metsIdx(i)));
    if not(isempty(formula))
        if strcmp(formula,'CO2') || strcmp(formula,'CO')
            CarbonAtoms(i,1) = 1;
        elseif strcmp(formula,'C2H')
            CarbonAtoms(i,1) = 2;
        else
            aux = formula(2:4);
            if not(isnan(str2double(aux)))
                CarbonAtoms(i,1) = str2double(aux);
            else
                aux = formula(2:3);
                if not(isnan(str2double(aux)))
                    CarbonAtoms(i,1) = str2double(aux);
                else
                    aux = formula(2:2);
                    if not(isnan(str2double(aux)))
                        CarbonAtoms(i,1) = str2double(aux);
                    else
                        CarbonAtoms(i,1) = 1;
                    end
                end
            end
        end
    else
        CarbonAtoms(i,1) = 0;
    end
end