% script to construct the ILP model, set an objective
% and generate the interactions based on this objective
clc
clear

modelDescription = {'myModel1','myModel2','myModel3','myModel4','...'};
% choose the indexes of the models we want to build an ILP model for
model_idx = [1,2];
%modelPath = 'yourPATH';
% saving directory should be the same us before
saving_directory = './output/';
%% making of the ILP model
model_ILP.species= modelDescription(model_idx);
model_ILP.DiMEsID = {};
model_ILP.DiMEs_secretions = {};
model_ILP.DiMEs_uptakes = {};
model_ILP.metabolites_all = {}; %all metabolites
model_ILP.metabolites = {}; %all metabolites of a certain species
model_ILP.varNames = {};
model_ILP.vartypes = {};
model_ILP.var_ub = [];
model_ILP.var_lb = [];
model_ILP.N=1;
model_ILP.M=1000;
model_ILP.description= concatenateList(modelDescription(model_idx),'|');


%for loop starts here to add variables

%Add DiMEsID under the format s#_l_#_alt_#

for index = 1:length(model_idx)
    
    table = load(strcat(saving_directory,modelDescription{model_idx(index)},'_summary.mat'));
    table = table.SummaryTable;
    load(strcat(saving_directory,modelDescription{model_idx(index)},'_unique_active.mat'),'tag');
    
    model_ILP.DiMEsID = [model_ILP.DiMEsID;table.Properties.RowNames];
    
    %Add secretions
    idx = find(contains(table.Properties.VariableNames,'P'));
    secretions = cell(size(table,1),1);
    
    for i=1:size(table,1)
        row = table(i,idx);
        add = find(~contains(row.Variables,{'none'}));
        met = row(1,add);
        met = met.Variables;
        secretions{i,1} = met;
    end
    
    model_ILP.DiMEs_secretions = [model_ILP.DiMEs_secretions;secretions];
    
    %Add uptakes
    idx = find(contains(table.Properties.VariableNames,'S'));
    uptakes = cell(size(table,1),1);
    for i=1:size(table,1)
        row = table(i,idx);
        add = find(~contains(row.Variables,{'none'}));
        met = row(1,add);
        met = met.Variables;
        uptakes{i,1} = met;
    end
    
    model_ILP.DiMEs_uptakes = [model_ILP.DiMEs_uptakes;uptakes];
    
    metabolites = cell(1,2);
    
    s=[secretions{:}];
    u=[uptakes{:}];
    
    metabolites{1,1} = unique(s);
    metabolites{1,2} = unique(u);
    
    model_ILP.metabolites = [model_ILP.metabolites ; metabolites];
    
    %Add all metabolites
    metabolites_all = table.Variables;
    metabolites_all = setdiff(metabolites_all,'none');
    
    model_ILP.metabolites_all = [model_ILP.metabolites_all;metabolites_all];
    model_ILP.metabolites_all = unique(model_ILP.metabolites_all,'stable');
    
    %ADD VARIABLES
    %Make y variables
    yields = unique(tag,'stable');
    for i=1:length(yields)
        model_ILP.varNames(end+1,1) = {strcat('y_s',num2str(index),'_',num2str(yields(i)))};
        model_ILP.var_ub(end+1,1) = 1;
        model_ILP.var_lb(end+1,1) = 0;
        model_ILP.vartypes(end+1,1) = {'B'};
    end
    
    %Make ws variables
    for i=1:size(metabolites{1,1},2)
        model_ILP.varNames(end+1,1) = strcat('ws','_s',num2str(index),'_',metabolites_all(strcmp(metabolites_all,metabolites{1,1}(1,i))));
        model_ILP.var_ub(end+1,1) = 1;
        model_ILP.var_lb(end+1,1) = 0;
        model_ILP.vartypes(end+1,1) = {'B'};
    end
    %Make wu variables
    for i=1:size(metabolites{1,2},2)
        model_ILP.varNames(end+1,1) = strcat('wu','_s',num2str(index),'_',metabolites_all(strcmp(metabolites_all,metabolites{1,2}(1,i))));
        model_ILP.var_ub(end+1,1) = 1;
        model_ILP.var_lb(end+1,1) = 0;
        model_ILP.vartypes(end+1,1) = {'B'};
    end
    
end

%Make z variables
    for i=1:size(model_ILP.DiMEsID,1)
        model_ILP.varNames(end+1,1) = strcat('z_', model_ILP.DiMEsID(i,1));
        model_ILP.var_ub(end+1,1) = 1;
        model_ILP.var_lb(end+1,1) = 0;
        model_ILP.vartypes(end+1,1) = {'B'};
    end

%Make xn and xp variables
for i=1:length(model_ILP.metabolites_all)
    model_ILP.varNames(end+1,1) = strcat('xn_', model_ILP.metabolites_all(i,1));
    model_ILP.var_ub(end+1,1) = 1;
    model_ILP.var_lb(end+1,1) = 0;
    model_ILP.vartypes(end+1,1) = {'B'};
end

for i=1:length(model_ILP.metabolites_all)
    model_ILP.varNames(end+1,1) = strcat('xp_', model_ILP.metabolites_all(i,1));
    model_ILP.var_ub(end+1,1) = 1;
    model_ILP.var_lb(end+1,1) = 0;
    model_ILP.vartypes(end+1,1) = {'B'};
end

%Make mu and nu variables
for i=1:length(model_ILP.metabolites_all)
    model_ILP.varNames(end+1,1) = strcat('mu_', model_ILP.metabolites_all(i,1));
    model_ILP.var_ub(end+1,1) = 1;
    model_ILP.var_lb(end+1,1) = 0;
    model_ILP.vartypes(end+1,1) = {'B'};
end

for i=1:length(model_ILP.metabolites_all)
    model_ILP.varNames(end+1,1) = strcat('nu_', model_ILP.metabolites_all(i,1));
    model_ILP.var_ub(end+1,1) = 1;
    model_ILP.var_lb(end+1,1) = 0;
    model_ILP.vartypes(end+1,1) = {'B'};
end

for i=1:length(model_ILP.metabolites_all)
    model_ILP.varNames(end+1,1) = strcat('u_', model_ILP.metabolites_all(i,1));
    model_ILP.var_ub(end+1,1) = 1;
    model_ILP.var_lb(end+1,1) = 0;
    model_ILP.vartypes(end+1,1) = {'B'};
end

model_ILP.constraintNames = {};
model_ILP.constraintType = {};
model_ILP.rhs = [];

%Objective function
model_ILP.f=zeros(size(model_ILP.varNames,1),1);
model_ILP.objtype = -1;

%Make A matrix
model_ILP.A=sparse(0,size(model_ILP.varNames,1));

for j=model_idx
    
    dimes = find(contains(model_ILP.DiMEsID,strcat(modelDescription{j},'_')));
    
    [~, NumVars] = size(model_ILP.A);
    
    %Add contraints on ws
    idx_z = find(contains(model_ILP.varNames,strcat('z_',model_ILP.species{j})));
    dime = strrep(model_ILP.varNames(idx_z),'z_','');
    dime = find(ismember(model_ILP.DiMEsID,dime));
    
    n=0;
    
    for i=1:length(idx_z)
        sec = model_ILP.DiMEs_secretions{dime(i),1};
        idx_ws = find(contains(model_ILP.varNames,strcat('ws_s',num2str(j))) & contains(model_ILP.varNames,sec));
        for k=1:length(idx_ws)
            NewCons = zeros(NumVars,1);
            NewCons(idx_z(i)) = 1;
            NewCons(idx_ws(k)) = -1;
            model_ILP.A(end+1,:) = NewCons;
            model_ILP.rhs(end+1,1) = 0;
            model_ILP.constraintNames{end+1,1} = strcat('SecAct_', num2str(n+k),'_s',num2str(j));
            model_ILP.constraintType{end+1,1} = '<';
        end
        n=n+k;
    end
    
    %Add contraints on wu
    n=0;
    
    for i=1:length(idx_z)
        up = model_ILP.DiMEs_uptakes{dime(i),1};
        idx_wu = find(contains(model_ILP.varNames,strcat('wu_s',num2str(j))) & contains(model_ILP.varNames,up));
        for k=1:length(idx_wu)
            NewCons = zeros(NumVars,1);
            NewCons(idx_z(i)) = 1;
            NewCons(idx_wu(k)) = -1;
            model_ILP.A(end+1,:) = NewCons;
            model_ILP.rhs(end+1,1) = 0;
            model_ILP.constraintNames{end+1,1} = strcat('UpAct_', num2str(n+k),'_s',num2str(j));
            model_ILP.constraintType{end+1,1} = '<';
        end
        n=n+k;
    end
    
    % Add constraints on sum of z and wu/ws
    
    [NumCons, ~] = size(model_ILP.A);
    idx_ws = find(contains(model_ILP.varNames,strcat('ws_s',num2str(j))));
    
    for i=1:length(idx_ws)
        %GET idx of DIMES WHERE THE METABOLITE IS PRESENT
        met = strrep(model_ILP.varNames(idx_ws(i)),strcat('ws_s',num2str(j),'_'),'');
        idx_dim = [];
        for k=dimes'
            if any(contains(model_ILP.DiMEs_secretions{k},met{1,1}))
                idx_dim = [idx_dim,k];
            end
        end
        idx_z_dim=find(ismember(model_ILP.varNames,strcat('z_',model_ILP.DiMEsID(idx_dim))));
        NewCons = zeros(NumVars,1);
        NewCons(idx_z_dim) = 1;
        NewCons(idx_ws(i)) = -1;
        model_ILP.A(NumCons+i,:) = NewCons;
        model_ILP.rhs(NumCons+i,1) = 0;
        model_ILP.constraintNames{NumCons+i,1} = strcat('Sec_', num2str(i),'_s',num2str(j));
        model_ILP.constraintType{NumCons+i,1} = '>';
    end
    %
    [NumCons, ~] = size(model_ILP.A);
    idx_wu = find(contains(model_ILP.varNames,strcat('wu_s',num2str(j))));
    
    for i=1:length(idx_wu)
        %GET idx of DIMES WHERE THE METABOLITE IS PRESENT
        met = strrep(model_ILP.varNames(idx_wu(i)),strcat('wu_s',num2str(j),'_'),'');
        idx_dim = [];
        for k=dimes'
            if any(contains(model_ILP.DiMEs_uptakes{k},met{1,1}))
                idx_dim = [idx_dim,k];
            end
        end
        idx_z_dim=find(ismember(model_ILP.varNames,strcat('z_',model_ILP.DiMEsID(idx_dim))));
        NewCons = zeros(NumVars,1);
        NewCons(idx_z_dim) = 1;
        NewCons(idx_wu(i)) = -1;
        model_ILP.A(NumCons+i,:) = NewCons;
        model_ILP.rhs(NumCons+i,1) = 0;
        model_ILP.constraintNames{NumCons+i,1} = strcat('Up_', num2str(i),'_s',num2str(j));
        model_ILP.constraintType{NumCons+i,1} = '>';
    end
    
    %Add constraints of y
    
    idx_y = find(contains(model_ILP.varNames,strcat('y_s',num2str(j))));
    
    [NumCons, ~] = size(model_ILP.A);
    
    for i=1:length(idx_y)
        yield_cut = strrep(model_ILP.varNames(idx_y(i)),strcat('y_s',num2str(j),'_'),'') ;
        idx_z_d = find(contains(model_ILP.varNames,strcat('z_',model_ILP.species(j),'_l_',num2str(yield_cut{1,1}),'_')));
        NewCons = zeros(NumVars,1);
        NewCons(idx_z_d) = 1;
        NewCons(idx_y(i)) = - model_ILP.N;
        model_ILP.A(NumCons+i,:) = NewCons;
        model_ILP.rhs(NumCons+i,1) = 0;
        model_ILP.constraintNames{NumCons+i,1} = strcat('yUse_', num2str(i),'_s',num2str(j));
        model_ILP.constraintType{NumCons+i,1} = '=';
    end
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_y) = 1;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = 1;
    model_ILP.constraintNames{NumCons+1,1} = strcat('yCut_',num2str(j));
    model_ILP.constraintType{NumCons+1,1} = '=';
    
end

%Add constraints that are for every j

for j=1:length(model_ILP.metabolites_all)
    
    met = model_ILP.metabolites_all{j};
    
    idx_xn = find(ismember(model_ILP.varNames,strcat('xn_',met)));
    idx_xp = find(ismember(model_ILP.varNames,strcat('xp_',met)));
    idx_mu = find(ismember(model_ILP.varNames,strcat('mu_',met)));
    idx_nu = find(ismember(model_ILP.varNames,strcat('nu_',met)));
    idx_wu = find(contains(model_ILP.varNames,'wu') & contains(model_ILP.varNames,strcat('_',met)));
    
    [NumCons, ~] = size(model_ILP.A);
    
    for i=1:length(idx_wu)
        idx_u = find(ismember(model_ILP.varNames,strcat('u_',met)));
        NewCons = zeros(NumVars,1);
        NewCons(idx_xp) = -1;
        NewCons(idx_u) = -1;
        NewCons(idx_wu(i)) = 1;
        model_ILP.A(NumCons+i,:) = NewCons;
        model_ILP.rhs(NumCons+i,1) = 0;
        model_ILP.constraintNames{NumCons+i,1} = strcat('Abiotic_', num2str(j),'_s',num2str(i));
        model_ILP.constraintType{NumCons+i,1} = '<';
    end
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_wu) = 1;
    NewCons(idx_xn) = -model_ILP.M;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = 1;
    model_ILP.constraintNames{NumCons+1,1} = strcat('UpComp1_',met);
    model_ILP.constraintType{NumCons+1,1} = '<';
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_wu) = -1;
    NewCons(idx_xn) = model_ILP.M;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = model_ILP.M - 2 ;
    model_ILP.constraintNames{NumCons+1,1} = strcat('UpComp2_',met);
    model_ILP.constraintType{NumCons+1,1} = '<';
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_wu) = 1;
    NewCons(idx_xp) = -1;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = 0;
    model_ILP.constraintNames{NumCons+1,1} = strcat('Coop1_',met);
    model_ILP.constraintType{NumCons+1,1} = '>';
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_wu) = 1;
    NewCons(idx_mu) = model_ILP.M;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = model_ILP.M;
    model_ILP.constraintNames{NumCons+1,1} = strcat('Coop2_',met);
    model_ILP.constraintType{NumCons+1,1} = '<';
    
    idx_ws = find(contains(model_ILP.varNames,'ws') & contains(model_ILP.varNames,strcat('_',met)));
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_ws) = 1;
    NewCons(idx_xp) = -1;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = 0;
    model_ILP.constraintNames{NumCons+1,1} = strcat('Coop3_',met);
    model_ILP.constraintType{NumCons+1,1} = '>';
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_ws) = 1;
    NewCons(idx_nu) = model_ILP.M ;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = model_ILP.M;
    model_ILP.constraintNames{NumCons+1,1} = strcat('Coop4_',met);
    model_ILP.constraintType{NumCons+1,1} = '<';
    
    [NumCons, ~] = size(model_ILP.A);
    
    NewCons = zeros(NumVars,1);
    NewCons(idx_mu) = 1;
    NewCons(idx_nu) = 1;
    NewCons(idx_xp) = 1;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1,1) = 1;
    model_ILP.constraintNames{NumCons+1,1} = strcat('Coop5_',met);
    model_ILP.constraintType{NumCons+1,1} = '>';
    
end

save(strcat(saving_directory,'model_ILP_',concatenateList(modelDescription(model_idx),'|'),'.mat'),'model_ILP')

%%
% setting an objective function
%[model_ILP_] = setObjective(model_ILP,'Cooperation','max');
[model_ILP_] = setObjective(model_ILP,'Competition','min');
% generate alternative interactions
maxAlt = 100;
solTime = 600;
filename = strcat(saving_directory,'ReMIND',concatenateList(modelDescription(model_idx),'|'),'.mat');
str = 'xn'; %variable used to generate the cut contsraints for the alternatives %xp ->cooperation, xn -> competition
[~,solutionTable] = ReMIND(model_ILP_,str,maxAlt,solTime,filename);
