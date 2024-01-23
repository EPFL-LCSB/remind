% .. Author:
% Anush Chiappino-Pepe 2018
% Evangelia Vayena 2023

%%
% FIRST RUN file_setup.m for each one pf the models
% and then run this script for all the models

%%%%%%%%%%%%%%%%%%%%%%
% DiMEs ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%
% inputs
if ~exist('grRate')
grRate = 0.2;%FBAsoln.f; % optimal growth rate (this can be obtained optimizing for growth)
end
essThr = 0.99; % essentiality threshold (% of optimal growth to be defined, see minObj). If a knockout leads to grow below this threshold the gene will be considered essential for growth.
minObj = essThr*grRate; %essThr*grRate; % minimal required growth
maxObj = grRate; % upper bound in growth
drainsForiMM = {}; % names of drains/exchanges used for IMM analysis. Empty means by default it will search for minimal media accross all drains/exchanges
metabData = []; % provide metabolomics data to integrate in this analysis. Empty if none.
NumAlt = 1000; % number of alternative minimal media to find. We define 1 for you to test that the pipeline works. But it is suggested to define 5000 or more to make sure you identify all alternatives
time = 600; % time limit for optimization in seconds. If empty we do not set a time limit.
tagMin = 0; % additional constrain to avoid generating suboptimal solutions, meaning that if 1 we will not identify media that is not minimal
tagBoth = 1; % 1 --> iMEs , 0 --> iMMs
tagOptYield = 1; % adding parsimonious uptake constraint
weight = 'CA'; %weights on parsimonious constraint 'MW' (molecular weight),'CA' (carbon atoms),'default'(all w=1)
tagCUT_pLOW=1; %add a lower bound on uptakes, this is used to constrain the solution space within a yield cut
tagMinSize=1; %if we want to enumerate all solutions of all sizes; 1 to add the block condition so that size of the alternatives <= min_size+2
filename = modeldescription;
p= [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
p= 1./p; % percentages for parsimonious, here you can discretise as you wish

% all carbon sources
[idx,carbonExchange] = getCarbonSources(model,'tfa'); %gets all carbon uptake reactions

for i = 1:length(carbonExchange)
    metsIdx(i,1) = find(model.S(:,carbonExchange(i)));
end

%Weights on parsiminious constraint
if strcmp(weight,'CA') % get carbon atoms
    w = getCarbonAtoms(model,metsIdx);
elseif strcmp(weight,'MW') % get MW(molecular weights)
    w = computeMW(model,model.mets(metsIdx));
elseif strcmp(weight,'default') %set weights as
    w = ones(length(idx),1);
end

%Constraint growth (f=objectve function)
model.var_ub(find(model.f)) = maxObj;
model.var_lb(find(model.f)) = minObj;

idxTFA = getCarbonSources(model,'tfa');
idx_C = find(ismember(model.rxns,strrep(model.varNames(idxTFA),'R_','')));

% Parsimonious uptake problem
model.f = zeros(length(model.varNames),1);
model.f(idx) = w; %set the weights for the exchange reactions

%Minimize
model.objtype = 1; % 1-> minmize , -1 -> maximize
sol_pUP = solveTFAmodelCplex(model,time); %solve to minimize flux of exchange reactions

model.f = zeros(length(model.varNames),1);
%Maximize
model.objtype = -1;
model.f(find(ismember(model.varNames,strcat('F_',model.rxns(find(model.c)))))) = 1; %put objective back to maximize growth

%if we do not want binary variables for the inorganics
f = find(ismember(model.rxns,inorganics));
idx_exchange2 = setdiff(idx_exchange,f);
%in case we want to keep track of some inorganics (e.g. nh4 if we are
%intersted in nitrogen sources and co2 secretion)
f = find(ismember(model.rxns,{'EX_nh4_e';'EX_co2_e'}));
idx_exchange2 = [idx_exchange2;f];


[model_p,drains,min_size] = analysisDiMEs(model,tagBoth,minObj, maxObj, model.rxns(idx_exchange2), ...
    metabData,time);

% Parsimonious uptake constraint
[NumCons, NumVars] = size(model_p.A);
NewCons = zeros(NumVars,1);
NewCons(idx) = w;
model_p.A(NumCons+1,:) = NewCons;
model_p.rhs(NumCons+1) = 1.01*p(1)*sol_pUP.val;
model_p.constraintNames{NumCons+1} = ['CUT_pUP'];
model_p.constraintType{NumCons+1} = '<';

if tagCUT_pLOW
    [NumCons, NumVars] = size(model_p.A);
    NewCons = zeros(NumVars,1);
    NewCons(idx) = w;
    model_p.A(NumCons+1,:) = NewCons;
    model_p.rhs(NumCons+1) = 0.99*p(1)*sol_pUP.val;
    model_p.constraintNames{NumCons+1} = ['CUT_pLOW'];
    model_p.constraintType{NumCons+1} = '>';
end

%here we can run sequentially for the different yield cuts but also in
%parallel (i.e. different matlab instances)
for i = 1:length(p)
    f = find(ismember(model_p.constraintNames,{'CUT_pUP'}));
    model_p.rhs(f) = 1.01*p(i)*sol_pUP.val;
    f = find(ismember(model_p.constraintNames,{'CUT_pLOW'}));
    model_p.rhs(f) = 0.99*p(i)*sol_pUP.val;
    
    
    %printConstraint(model_p,{'CUT_pUP'})
    fprintf('Generating alternatives for %4.2f % of the optimal biomass yield \n',1/p(i)*100)
    saving_directory_new = strcat(strcat(saving_directory,num2str(1/p(i)*100),'%_opt_yield/'));
    mkdir(saving_directory_new);

    sol = solveTFAmodelCplex(model_p);
    if ~isempty(sol.val)
        
        [DPsimm,modelmilp2] = findDPMax_DiMEs(model_p, NumAlt, model_p.indUSE,time, tagMin, strcat(saving_directory_new,filename),tagMinSize);
        
        % save workspace
        save(strcat(saving_directory_new,filename,'.mat'));
        
        % save results
        filename = modeldescription;
        r1 = load(strcat(saving_directory_new,filename,'.mat'));
        
        if not(isempty(DPsimm))
            % Extract info about composition of the IMMs
            [immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
                immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfo_DiMEs(r1.model_p, r1.DPsimm, ...
                r1.model_p.indUSE,  r1.model_p.rowDrain_F,...
                r1.model_p.rowDrain_R,r1.model_p.F_Flux,  r1.model_p.R_Flux) ;
            
            
            immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
                num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];
            
            writeData(strcat(saving_directory_new,filename,'_',num2str(1/p(i)*100),'.csv'), immOutput.summary,...
                '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
                'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
                '%s\t%s\t%s\t%s\t%s');
            
            save(strcat(saving_directory_new,'Active_',modeldescription,'.mat'),'Active');
            clear Active
            
            elimIntermPMFile(saving_directory_new,filename)
        end
        f = find(contains(model_p.constraintNames,'CUTimm'));
        model_p.constraintNames(f) = strcat(model_p.constraintNames(f),'_old');
    end
end

% Get unique solutions and prepare the input for the ILP step of ReMIND
getUniqueSol(p,modeldescription,saving_directory)
run GenerateSummaryTables.m