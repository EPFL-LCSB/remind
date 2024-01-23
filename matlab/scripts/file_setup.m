%% Setting up the paths, solver and preparing the model for ReMIND
clear
clc

% Provide path to folder where your results will be saved
saving_directory = 'MyPath/output/';
if ~isdir(saving_directory)
    mkdir(saving_directory)
end

%add paths
ReMIND_directory =  'MyPaths/matReMIND/';
mattfa_directory = 'MyPath/mattfa/';
cplex_directory =  'MyPathToCplex';
addpath(genpath(ReMIND_directory));
addpath(genpath(mattfa_directory));
addpath(genpath(cplex_directory));
changeCobraSolver('cplex_direct');

% Tag true to apply phenomapping with TFA (thermodynamic constraints will
% be taken into account)
tagThermo = 0;

% Prepare model for phenomapping
%load model and thermo data

modelPath = 'MyPathToModel/modelName.mat';
modeldescription = strrep(modelPath,'MyPathToModel/','');
modeldescription = strrep(modeldescription,'.mat','');
load(modelPath)

%% set envirnoment
%Define the envirnoment you want to consider for the ReMIND workflow
% block uptakes froom demand and exchange reactions
% allow only the uptakes you want to consider for the analysis
idx_exchange = find(findExcRxns(model));
model.lb(idx_exchange) = 0;

%for example
inorganics = {'EX_cl_e';'EX_ca2_e';'EX_cobalt2_e';'EX_mobd_e';...
    'EX_cu2_e';'EX_fe2_e';'EX_fe3_e';'EX_h2o_e';'EX_k_e';'EX_mg2_e';...
    'EX_mn2_e';'EX_na1_e';'EX_pi_e';'EX_so4_e';...
    'EX_zn2_e';'EX_o2_e';'EX_h_e';'EX_nh4_e'};
f = find(ismember(model.rxns,inorganics));
model.lb(f) = -25;


% vitamins and carbon source
vitamins = {'EX_thm_e';'EX_btn_e';...
    'EX_ribflv_e';'EX_pydx_e';'EX_4abz_e';'EX_lipoate_e';...
    'EX_nac_e';'EX_fol_e';'EX_cbl1_e';'EX_pnto__R_e'};
f = find(ismember(model.rxns,vitamins));
model.lb(f) = -1;

% here we can define the main carbon sources
carbonSources = {'EX_fru_e';'EX_glc__D_e';'EX_ac_e'};
f = find(ismember(model.rxns,carbonSources));
model.lb(f) = -10;

% Bound secretions to get realistic DiMEs
% To bound all secretions
model.ub(idx_exchange) = 25;

% block sinks
f = find(contains(model.rxns,'sink_'));
model.ub(f) = 0;
% block the secretion of certain metabolites if necessary
f = find(ismember(model.rxns,{'EX_glc__D_e'}));
model.ub(f) = 0;
% if growth rate if known
growthRate = 0.2;

%%
if tagThermo ==  0
    [model, checkList, tagReady] = initDiMEsModel(model,'',tagThermo);
else
    % this is the thermo database included in the mattfa toolbox
    thermo_data_directory = '/yourPATH/mattfa/thermoDatabases/thermo_data.mat';
    thermo_data = load(thermo_data_directory);
    thermo_data = thermo_data.DB_AlbertyUpdate;
    [model, checkList, tagReady] = initDiMEsModel(model,thermo_data,tagThermo);
end
% Save prepared model
save(strcat(saving_directory,modeldescription,'_4ReMIND.mat'), 'model');