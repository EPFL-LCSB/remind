function [active,solutionTable] = ReMIND(model_ILP,str,maxAlt,solTime,filename)

%[model_ILP,str] = setObjective(model_ILP,obj,type);

Sol = solveTFAmodelCplex(model_ILP,solTime);

if (Sol.val==0)
    fprintf('Sol.val = 0\n')
end

active = []; %full solution

Alternative = [];
Length=[];
Metabolites = {};
Uptake = {};
Secretion = {};
Dimes = {};

numAlt=1;

while (~(isempty(Sol.x)) && (numAlt <= maxAlt)) 

    fprintf('Alternative : %i\n',numAlt)

    active(:,numAlt) = Sol.x;
    Alternative(numAlt,1) = numAlt;
    Length(numAlt,1) =  Sol.val;

    x = getAllVar(model_ILP,{str});
    idx = x(Sol.x(x)==1);
    ac = model_ILP.varNames(idx)';

    Metabolites{numAlt,1} = strrep(ac,strcat(str,'_'),'');

    xp = getAllVar(model_ILP,{'xp'});
    idxp = xp(Sol.x(xp)==1);
    idxpi= xp(Sol.x(xp)==0);
    xn = getAllVar(model_ILP,{'xn'});
    idxn = xn(Sol.x(xn)==1);
    idxni= xn(Sol.x(xn)==0);
    z = getAllVar(model_ILP,{'z'});
    idxz = z(Sol.x(z)==1);

    [NumCons, NumVars] = size(model_ILP.A);

    NewCons = zeros(NumVars,1);
    NewCons(idxz) = 1;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1) = length(idxz)-1;
    model_ILP.constraintNames{NumCons+1} = strcat('Alt_cut_',num2str(numAlt));
    model_ILP.constraintType{NumCons+1} = '<';

    [NumCons, NumVars] = size(model_ILP.A);

    NewCons = zeros(NumVars,1);
    NewCons(idxp) = -1;
    NewCons(idxn) = -1;
    NewCons(idxpi) = 1;
    NewCons(idxni) = 1;
    model_ILP.A(NumCons+1,:) = NewCons;
    model_ILP.rhs(NumCons+1) = 1-length(idxp)-length(idxn);
    model_ILP.constraintNames{NumCons+1} = strcat('Alt_',num2str(numAlt));
    model_ILP.constraintType{NumCons+1} = '>';
    
    wu = getAllVar(model_ILP,{'wu'});
    idx = wu(Sol.x(wu)==1);
    wus = model_ILP.varNames(idx);
    up = wus(find(contains(wus,strcat('_',Metabolites{numAlt,1}))));
    Uptake(numAlt,1) = {up'};

    ws = getAllVar(model_ILP,{'ws'});
    idx = ws(Sol.x(ws)==1);
    wus = model_ILP.varNames(idx);
    up = wus(find(contains(wus,strcat('_',Metabolites{numAlt,1}))));
    Secretion(numAlt,1) = {up'};

    dim = model_ILP.varNames(idxz)';
    Dimes(numAlt,1) = {dim'};

    Sol = solveTFAmodelCplex(model_ILP,solTime);
    numAlt=numAlt+1;
end

solutionTable = table(Alternative,Length,Metabolites,Uptake,Secretion,Dimes);
save(filename,'active','solutionTable')