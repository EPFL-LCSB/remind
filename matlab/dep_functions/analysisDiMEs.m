function [model, drains,min_size] = analysisDiMEs(model,exchange, minObj, ...
    maxObj, drainsForiMM, metabData,time)
% Identifies In silico Minimal Media (IMM) or In silico Minimal Secretion
% (IMS) to achieve a given objective
%
% USAGE:
%
%    [model, drains] = analysisIMM(model, flagUpt, minObj, maxObj, drainsForiMM, rxnNoThermo, ReactionDB, metabData)
%
% INPUT:
%    model:           TFA model structure
%
% OPTIONAL INPUTS:
%    flagUpt:         True to identify the In silico Minimal Media (IMM).
%                     Else it gets the In silico Minimal Secretion (IMS).
%                     (default = true)
%    minObj:          Objective lower bound (default = 90% optimal value)
%    maxObj:          Objective upper bound (default = optimal value)
%    drainsForiMM:    Drains or transports to minimize (default = all drains)
%    metabData:       Metabolomics data (default = empty)
%
%
% OUTPUTS:
%    model:           Model with MILP formulation for IMM analysis
%    drains:          Drains used for the IMM analysis
%    modelpre:        Model ready for MILP formulation
%
% .. Author:
% Meric Ataman 2014 : initial problem formulation
% Anush Chiappino-Pepe 2017 : minimal secretion and refinement of function
%

solWT = solveTFAmodelCplex(model,time);


if (nargin < 2)
    minObj = 0.9*solWT.val;
end
if (nargin < 3)
    maxObj = solWT.val;
end
if (nargin < 4)
    drainsForiMM = {};
end
if (nargin < 5)
    metabData = [];
end

if (minObj > maxObj)
    error('the defined lower bound (minObj) is bigger than the upper bound (maxObj) of the objective')
end

fprintf('getting model drains\n');
[model, flagChange, drains, drainMets] = putDrainsForward(model);
drainMets = {};
aux = findExcRxns(model);
drains = model.rxns(find(aux));
for i = 1:length(drains)
    f = find(ismember(model.rxns,drains{i}));
    drainMets(i,1) = model.mets(find(model.S(:,f)));
end
if ~isempty(drainsForiMM)
    if ((sum(ismember(drainsForiMM,model.rxns)) == length(drainsForiMM)) || (sum(ismember(drainsForiMM,drains)) == length(drainsForiMM)))
        drains = drainsForiMM;
        drainMets = printRxnFormula(model,drains,0,0,1);
    else
        fprintf('CAUTION: Not all drainsForiMM were identified as drains or rxns\n');
        fprintf('The analysis will be done for all substrates\n');
    end
end


drains_F = [ strcat('F_', drains) drainMets]; % get exchange
drains_R = [strcat('R_', drains)  drainMets]; % get exchange


if flagChange
    fprintf('some drains need to be redefined -> reconvert to thermo\n');
    error('please run initTestPhenoMappingModel before calling this function')
end
if ~isempty(metabData)
    model = loadConstraints(model,metabData);
end

model.var_lb(model.f==1) = minObj;
model.var_ub(model.f==1) = maxObj;
model.f = zeros(length(model.varNames),1);
model.objtype = -1; %maximize

% constrain binary variableS

[~,rowDrain_F ] = ismember(strcat('FU_',drains(:,1)), model.varNames);
[~,rowDrain_R ] = ismember(strcat('BU_',drains(:,1)), model.varNames);
% % % constrain continues variableS
%  [~,rowDrain_F ] = ismember(strcat('F_',drains(:,1)), model.varNames);
%  [~,rowDrain_R ] = ismember(strcat('R_',drains(:,1)), model.varNames);

fprintf('defining MILP problem\n');
intTag = {'BFUSE'}; % integer tag used in this MILP
[~,num_vars] = size(model.A);
indUSE = zeros(size(drains,1),1);
for i=1:size(drains,1)
    model.varNames(num_vars+i) = strcat(intTag, '_', drains(i,1));
    model.var_ub(num_vars+i) = 1;
    model.var_lb(num_vars+i) = 0;
    model.vartypes(num_vars+i) = {'B'};
    model.f(num_vars+i) = 1;
    indUSE(i) = num_vars+i;
end
model.indUSE = indUSE;
[~,model.rowDrain_F ] = ismember(strcat('FU_',drains(:,1)), model.varNames);
[~,model.rowDrain_R ] = ismember(strcat('BU_',drains(:,1)), model.varNames);

[~,model.F_Flux ] = ismember(drains_F(:,1), model.varNames);
[~,model.R_Flux ] = ismember(drains_R(:,1), model.varNames);

% define constraints for MILP
% 50 * BFUSE_R_rxn + R_rxn + F_rxn < 50
% BFUSE_R_rxn = 1 -> R_rxn + F_rxn < 0
% BFUSE_R_rxn = 0 -> R_rxn + F_rxn < 50

% BFUSE_R_rxn + R_rxn + F_rxn > 10^-7
% BFUSE_R_rxn = 1 -> R_rxn + F_rxn > 0
% BFUSE_R_rxn = 0 -> R_rxn + F_rxn > 10^-7

% constrain binary variableS
if exchange == 0
    [num_constr,~] = size(model.A);
    for i=1:length(rowDrain_R)
        model.rhs(num_constr+i,1) =  1000;
        model.constraintNames{num_constr+i,1} = strcat('BFMER_',num2str(i));
        model.constraintType{num_constr+i,1} = '<';
        model.A(num_constr+i,rowDrain_R(i)) = 1;
        model.A(num_constr+i,num_vars+i) = 1000;
    end
    [num_constr,~] = size(model.A);
    
    for i=1:length(rowDrain_R)
        model.rhs(num_constr+i,1) =  1e-07;
        model.constraintNames{num_constr+i,1} = strcat('BFMER1_',num2str(i));
        model.constraintType{num_constr+i,1} = '>';
        model.A(num_constr+i,rowDrain_R(i)) = 1;
        model.A(num_constr+i,num_vars+i) = 1;
    end
else
    [num_constr,~] = size(model.A);
    for i=1:length(rowDrain_R)
        model.rhs(num_constr+i,1) =  1000;
        model.constraintNames{num_constr+i,1} = strcat('BFMER_',num2str(i));
        model.constraintType{num_constr+i,1} = '<';
        model.A(num_constr+i,rowDrain_F(i)) = 1;
        model.A(num_constr+i,rowDrain_R(i)) = 1;
        model.A(num_constr+i,num_vars+i) = 1000;
    end
    [num_constr,~] = size(model.A);
    
    for i=1:length(rowDrain_R)
        model.rhs(num_constr+i,1) =  1e-05;
        model.constraintNames{num_constr+i,1} = strcat('BFMER1_',num2str(i));
        model.constraintType{num_constr+i,1} = '>';
        model.A(num_constr+i,rowDrain_F(i)) = 1;
        model.A(num_constr+i,rowDrain_R(i)) = 1;
        model.A(num_constr+i,num_vars+i) = 1;
    end
end

sol = solveTFAmodelCplex(model,time);
min_size = length(model.indUSE) - sol.val;
end

