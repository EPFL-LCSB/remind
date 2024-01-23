function [f,carbonExchange] = getCarbonSources(model,mode)    
%function that get all carbon uptake reactions

ExtraFormulas = {'CNO','CHO3','CNS','CHN'};
NoFormulaStr = {'NA','None','NULL',''};
InorganicMetID = find(cellfun(@isempty,regexp(model.metFormulas,'(?-i)C[A-Z_0-9]'))); %finds inorg metabolites with carbons
NoFormulaStrID = find(ismember(model.metFormulas,NoFormulaStr)); %find metabolites with no formula id
%InorganicMetID = setdiff(InorganicMetID,NoFormulaStrID);
InorganicMetID = unique([InorganicMetID;NoFormulaStrID]); %removes repetitions
ExtraFormulasID = find(ismember(model.metFormulas,ExtraFormulas)); %find the extra formulas met id
InorganicMetID = union(InorganicMetID,ExtraFormulasID); %combines the two

%Find exchange reactions in the model
f1 = findExcRxns(model);
f1 = find(f1);

for i = 1:length(f1) %for all exchange reactions
            metsIdx = find(model.S(:,f1(i))); %find the metabolites id involved in the exch react 
            %from stch matrix (we know the reaction (columns), so we get the non-zero indx in all rows(metabolites))
            if strcmp(model.metCompSymbol(metsIdx),'c') %find metabolite with c localization (in the cell)
                f1(i) = 0;
            end
        end
f1 = f1(find(f1)); %consider only exch reaction not-involving metabolites with c localization
%exclude sinks (c is for contained, so it removes the intracellular metabolites)
d = 0;

for i = 1:length(f1) 
    f2 = find(model.S(:,f1(i))); %find the metabolites id involved in the exch react (non-c)
    if ismember(f2,InorganicMetID) 
    else %if the met is not in the list of inorganic met
        d = d + 1;
        f3(d,1) = f1(i); %add to f3
    end
end
carbonExchange = f3;
if strcmp(mode,'tfa') %if the mode was set as tfa
   f3 = strcat('R_',model.rxns(f3)); %add R_ to the rxn to indicate uptake
   f = find(ismember(model.varNames,f3)); %find the R reactions
else
    f = f3; 
end

end
