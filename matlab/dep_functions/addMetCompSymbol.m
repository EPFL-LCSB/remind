function model = addMetCompSymbol(model)
% Adds the field metCompSymbol to a model
%
% USAGE:
%
%    model = addMetCompSymbol(model)
%
% INPUT:
%    model:           FBA/TFA model structure
%
% OUTPUTS:
%    model:           Model with field metCompSymbol
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% check if the compartments for the metabolites are defined as in the RAVEN
% toolbox: the conversion is then easy
if isfield(model,'comps') && isfield(model,'metComps') && ...
        iscell(model.comps) && isfloat(model.metComps)
    model.metCompSymbol = model.comps(model.metComps);
else
    % we consider the following descriptions to assign compartments to mets:
    % _compartmentID or [compartmentID] or (compartmentID) at the end of the
    % mets field in the model
    symbols = {'_','[','('};
    
    % get three last characters of model.mets field
    metCompTag = cell(length(model.mets),1);
    for i = 1:length(model.mets)
        metCompTag{i} = model.mets{i}(end-2:end);
    end
    
    % check what tags are part of the three last characters
    numTag = zeros(length(model.mets),length(symbols));
    for i = 1:length(symbols)
        numTag(:,i) = ~cellfun(@isempty,regexpi(metCompTag,symbols{i},'match'));
    end
    
    % assign compartment to metabolite based on tag
    metCompSymbol = cell(length(model.mets),1);
    for i = 1:length(model.mets)
        if sum(numTag(i,:)) == 0
            warning('one met does not have a compartment tag. It was assigned to the cytosol. Please check the metabolite, the compartment tags of this function and verify this assumption is correct')
            fprintf(strcat('Problematic metabolite ID:',model.mets{i},'\n'));
            metCompSymbol{i} = 'c';
        elseif sum(numTag(i,:)) > 1.5
            warning('one met has more than one compartment tag at the end. WTF? We assinged the last letter of the string as the compartment. Please check the metabolite and the compartment tags of this function and verify this assumption is correct')
            fprintf(strcat('Problematic metabolite ID:',model.mets{i},'\n'));
            metCompSymbol{i} = metCompTag{i}(end);
        elseif sum(numTag(i,:)) < 1.5
            if strcmp(symbols(numTag(i,:)>0.5),'_')
                metCompSymbol{i} = metCompTag{i}(end);
            else
                metCompSymbol{i} = metCompTag{i}(end-1);
            end
        end
    end
    model.metCompSymbol = metCompSymbol;
end
