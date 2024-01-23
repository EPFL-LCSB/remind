function [model, checkList, tagReady] = initDiMEsModel(modeli, ...
    ReactionDB, tagThermo, rxnNoThermo)
% Initial check list to verify the model is ready for a PhenoMapping
% analysis
%
% USAGE:
%
%       [model, checkList, tagReady] = initTestPhenoMappingModel(modeli, ReactionDBpath, tagThermo, rxnNoThermo)
%
% INPUTS:
%    model:           model with FBA/TFA structure
%
% OPTIONAL INPUTS:
%    ReactionDBpath:  Path to database with all the thermodynamic
%                     information from Alberty (default = empty / to
%                     provide manually)
%    tagThermo:       True if it is desired to convert the model to thermo
%                     with the matTFA toolbox (default = true)
%    rxnNoThermo:     List of rxns for which no thermo constraints should
%                     be applied (default = empty / no reaction relaxed)
%
%
% OUTPUTS:
%    model:           model ready for PhenoMapping, unless it is missing
%                     major pieces of information and it requires a manual
%                     curation
%    checkList:       list of tests done
%    tagReady:        True if model is ready for PhenoMapping
%
% Anush Chiappino-Pepe 2018
%Evangelia Vayena 2022

if (nargin < 2)
    ReactionDB = [];
end
if (nargin < 3)
    tagThermo = 1;
end
if (nargin < 4)
    rxnNoThermo = [];
end

tagReady = 1;
model = modeli;
checkList = cell(5,1);

% test for feasibility of the model
sol = solveFBAmodelCplex(model);
if isnan(sol.f) || isempty(sol.f) || sol.f<1E-3
    error('the model provided is not feasible!')
end

% check the model structure required for thermo conversion
fprintf('1: checking model structure\n');
if ~isfield(model,'metCompSymbol')
    warning('the model does not contain the field metCompSymbol. The compartments will be searched automatically')
    model = addMetCompSymbol(model);
end

if tagThermo && (~isfield(model,'metSEEDID') || ~isfield(model,'CompartmentData'))
    warning('the model does not contain the fields metSEEDID or/and CompartmentData required in matTFA for conversion to thermo. Hence it wont be converted to thermo')
    tagThermo = 0;
end

% check that drains are all of the type A => (required for in silico
% minimal media analysis)
[model, flagChange] = putDrainsForward(model);
if flagChange
    checkList{1} = 'corrected: some drains in the model were not defined as A=> and the model was reconverted to a TFA structure';
else
    checkList{1} = 'ok: the model had all drains of the type A=>';
end

% check that the model has a TFA-friendly structure
if isfield(model,'A') && isfield(model,'f') && isfield(model,'var_lb') ...
        && isfield(model,'var_ub') && isfield(model,'rhs') ...
        && isfield(model,'constraintNames') && isfield(model,'varNames') ...
        && isfield(model,'constraintType') && isfield(model,'vartypes') ...
        && isfield(model,'objtype')
    checkList{2} = 'ok: the model already has a TFA structure';
    tagConv = 0;
else
    tagConv = 1;
    checkList{2} = 'corrected: the model was converted to a TFA structure';
end

% if the model does not have a TFA-friendly structure or it has drains like
% =>A create a TFA structure
if flagChange || tagConv
    if tagThermo
        
        modelt = prepModelforTFA(model, ReactionDB, model.CompartmentData);
        model = convToTFA(modelt, ReactionDB, [], 'DGo', {}, 0.001, 1, 1);
        
        % it can happen that a model is not thermodynamically feasible because of
        % the directionality of a set of reactions. One should investigate these
        % cases in more detail. One option is to remove the thermodynamic
        % constraints for these reactions: this will be the input "rxnNoThermo".
        % Here, if the model is thermo infeasible, we will only apply FBA.
        sol = solveTFAmodelCplex(model);
        if isnan(sol.val) || isempty(sol.val) || sol.val<1E-3
            warning('the model is not thermodynamically feasible: we took out the thermo constraints. To work with thermo, please, investigate the reactions that make your model infeasible with the matTFA repository')
            model = convGEMToTFAStruc(model);
        end
    else
        model = convGEMToTFAStruc(model);
    end
end


sol = solveTFAmodelCplex(model);
if isnan(sol.val) || isempty(sol.val) || sol.val<1E-3
    error('the model is not feasible')
else
    fprintf('the model is feasible!\n');
end

% get the net fluxes associated with TFA structure
model.indNF = getAllVar(model,{'NF'});
if isempty(model.indNF)
    model = addNetFluxVariables(model);
    model.indNF = getAllVar(model,{'NF'});
end

% make sure the model has the structures required for gene essentiality
fprintf('2: checking fields for gene essentiality\n');
if isfield(model,'genes') && isfield(model,'grRules')
    checkList{3} = 'ok: the model has the fields for gene essentiality';
    if isfield(model,'rules')
        checkList{4} = 'ok: the field rules is present';
    else
        [model] = generateRules(model);
        checkList{4} = 'corrected: the field rules was added';
    end
else
    checkList{3} = 'issue: add fields genes and grRules to the model';
    checkList{4} = 'not tested: presence of field rules';
    tagReady = 0;
end

% summarize final status of the model
if isequal(rmfield(model,'indNF'),modeli) && tagReady
    fprintf('all checks passed: the output model is the same as the input model. The model is ready for phenomapping\n');
elseif ~isequal(rmfield(model,'indNF'),modeli) && tagReady
    fprintf('checks were applied: the output model is now ready for phenomapping\n');
else
    fprintf('some checks failed: the output model needs manual curation for phenomapping - details in checkList\n');
end