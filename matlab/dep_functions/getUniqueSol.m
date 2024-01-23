function [] = getUniqueSol(p,modelname,results_dir)

%GEt unique DiMEs keeping only highest yields

%Name of the file to open
file_name = strcat('/',modelname);

A = [];
tag = [];

for i = 1:length(p)
    if exist(strcat(results_dir,num2str(1/p(i)*100),'%_opt_yield',file_name,'_',num2str(1/p(i)*100),'.csv'))
       csv_table = readtable(strcat(results_dir,num2str(1/p(i)*100),'%_opt_yield',file_name,'_',num2str(1/p(i)*100),'.csv'));
       Drains = strrep(csv_table.DrainID,'EX_','');
       break
    end
end

for i = 1:length(p)
    if exist(strcat(results_dir,num2str(1/p(i)*100),'%_opt_yield','/Active_',modelname,'.mat'))
        load(strcat(results_dir,num2str(1/p(i)*100),'%_opt_yield','/Active_',modelname,'.mat'));
        A = [A Active];
        tag = [tag ; (1/p(i)*100).*ones(size(Active,2),1)];
    end
end

    [Active,idx] = unique(A','rows','stable');  
    Active = Active';
    tag = tag(idx);
    save(strcat(results_dir,modelname,'_unique_active.mat'),'Active','tag','Drains')
end

