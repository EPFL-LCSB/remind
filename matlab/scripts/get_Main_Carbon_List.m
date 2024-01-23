%this script is used to differantiate between main carbon sources and
%compounds that are neceessary due to auxtrophies
%More specifically, we identify metabolites that can act as the main carbon
%source (big uptake flux) and others that are used to deal with
%auxotrophies (small uptake flux)
%we can already have a list of compounds that we want to consider as main
%carbon sources but maybe we want to expand it or
%for some models, the list with the main carbon sources we have may not be
%sufficient, thus this script will also enrich our main carbon list.
%%
clc
clear

%add paths
ReMIND_directory =  '/yourPATH';
mattfa_directory = '/yourPATH';
cplex_directory =  '/yourPATH';
addpath(genpath(ReMIND_directory));
addpath(genpath('/'));

saving_directory = '/yourPATH';
if ~isdir(saving_directory)
    mkdir(saving_directory)
end

%load model and thermo data
modeldescription = 'myModel';
modelPath = ['yourPATH/matReMIND/models/myModel.mat'];
load(modelPath)

% here the set up is for the tfa model
% in case you don't want to use thermodynamic constraints
% you can use the line below to convert it in the tfa structure
model = convGEMToTFAStruc(model);

%Add metCompSymbol field - needed for getCarbonSources.m
if ~isfield(model,'metCompSymbol')
    f1 =find(contains(model.mets,'_c'));
    model.metCompSymbol(f1,1) = {'c'};
    f1 =find(contains(model.mets,'_e'));
    model.metCompSymbol(f1,1) = {'e'};
    f1 =find(contains(model.mets,'_p'));
    model.metCompSymbol(f1,1) = {'p'};
end


% INPUTS
main = {'EX_glc__D_e';'EX_fru_e'};% this is the list of the main carbon sources we consider for the analysis
rho = 3;        %small flux coeff
M = 25;         %big flux coeff
TagEss = 0;     %if this is 1, the essential carbons sources will be blocked
grRate = 0.2;   % optimal growth rate (this can be obtained optimizing for growth)
essThr = 0.99;  % essentiality threshold (% of optimal growth to be defined, see minObj).
minObj = essThr*grRate;   % minimal required growth
maxObj = grRate;          % upper bound in growth
MaxAlt = 100;    %max # of alternatives to look for
tagAerobic = 0; %1 for aerobic, 0 for anaerobic conditions
tagExpand = 1; %1 if we just want to expand the current carbon list
N = 1; %  parameter for catabolite repression - how many main carbon sources we allow to be uptaken simultaneously
filename = strcat(modeldescription,'_main_carbon_list');
pathsave = saving_directory;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Find all Carbon Sources
idxTFA = getCarbonSources(model,'tfa');
idx_C = find(ismember(model.rxns,strrep(model.varNames(idxTFA),'R_','')));

%Find C sources from the Main list present in the model
main = main(find(ismember(main,model.rxns)));
%the remaining carbon sources in the model can be used to support growth
Rest = setdiff(model.rxns(idx_C),main);

idx = find(findExcRxns(model));
idx_sink = find(contains(model.rxns,{'sink_';'DM_'}));
idx = setdiff(idx,idx_sink);

%Constraint growth
model.var_ub(find(model.f)) = maxObj;
model.var_lb(find(model.f)) = minObj;

%Set bounds on each reaction
idx_R = find(ismember(model.varNames,strcat('R_',model.rxns(idx))));
idx_F = find(ismember(model.varNames,strcat('F_',model.rxns(idx))));

model.var_ub(idx_F) = 25;
model.var_ub(idx_R) = 25;

%So it doesn't secrete glucose or fructose
if ismember('EX_glc__D_e',model.rxns)
    f=find(ismember(model.varNames,'F_EX_glc__D_e'));
    model.var_ub(f) = 0;
end

if ismember('EX_fru_e',model.rxns)
    f=find(ismember(model.varNames,'F_EX_fru_e'));
    model.var_ub(f) = 0;
end

idx_Rr = find(ismember(model.varNames,strcat('R_',Rest)));
idx_Rm = find(ismember(model.varNames,strcat('R_',main)));

model.var_ub(idx_Rm) = M;
model.var_ub(idx_Rr) = rho;

if tagAerobic
    o2_idx=find(ismember(model.varNames,'R_EX_o2_e'));
    model.var_ub(o2_idx)=20;
else
    o2_idx=find(ismember(model.varNames,'R_EX_o2_e'));
    model.var_ub(o2_idx)=0;
end

model = catabolite_repression(model,main,N);
% check if the model is feasible with these settings
sol = solveTFAmodelCplex(model);
%%
% if not, find the minimum number of metabolites to add to the main carbon
% list - here we can have some metabolites that we don't want to add to the
% main carbon list (e.g. they are essential for growth or irrelevant to our
% study ( see bellow Ess_list
if isempty(sol.val) || isnan(sol.val) || tagExpand
    %Add variables for the MILP problem
    %MFA
    model = addMFAvar(model,Rest);
    %BFUSE
    model = addBFUSE(model,Rest);
    
    %Add main constraints
    model = constrain_uptake(model,Rest,rho,M);
    
    %Block or not essentials
    if TagEss
        Ess_list = {'EX_uri_e';'EX_cytd_e'};
        model = block_essential(model,Ess_list);
    end
    
    % Get indexes
    [~,idx_BFUSE] = ismember(strcat('BFUSE_',Rest),model.varNames);
    [~,idx_BU] = ismember(strcat('BU_',Rest),model.varNames);
    [~,idx_MFA] = ismember(strcat('MFA_',Rest),model.varNames);
    
    %Define objective function as maximizing the sum of BFUSE ( so
    %minimizing the uptakes)
    
    [~,NumVars] = size(model.A);
    
    model.f = zeros(NumVars,1);
    model.f(idx_BFUSE) = 1;
    model.objtype = -1;
    
    %Solve
    model.var_ub(idx_Rr) = M;
    sol = solveTFAmodelCplex(model);
    
    if ~isempty(sol.val)
    %Get alternatives
    [model,main_c,DPs]=find_altern(model,MaxAlt,idx_BFUSE,idx_BU,idx_MFA);
    end % here if we dont find a solution we can increase the N and/or change the values of rho and  M
    
    main = [main;main_c];
    
    T=cell2table(main);
    writetable(T,strcat(saving_directory,filename,'.csv'));
    
end
