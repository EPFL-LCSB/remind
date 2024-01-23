%Define the parsimonious values we used
rowsP=[];
rowsS=[];
rownames={};

load(strcat(saving_directory,modeldescription,'_unique_active.mat'))
Drains = Drains';
for i = 1:size(Active,2)
max_P(i) = length(find(Active(:,i)==1));
end
for i = 1:size(Active,2)
max_S(i) = length(find(Active(:,i)==-1));
end
max_P = max(max_P);
max_S= max(max_S);
P = cell(max_P,1)';
S = cell(max_S,1)';

for j=1:size(Active,2)
    idx_p_j=find(Active(:,j)==1);
    idx_s_j=find(Active(:,j)==-1);
    newrowsP=[Drains(idx_p_j)];
    newrowsS=[Drains(idx_s_j)];
    P(j,1:length(idx_p_j)) = newrowsP;
    S (j,1:length(idx_s_j)) = newrowsS;
    rownames= horzcat(rownames,strcat(modeldescription,'_l_',num2str(tag(j)),'_alt_',num2str(j)));
    newrowsP={};
    newrowsS={};
end


columnNamesP={};
columnNamesS={};
for i= 1:size(P,2)
    columnNamesP{i}=strcat('P_',int2str(i));
end
for i= 1:size(S,2)
    columnNamesS{i}=strcat('S_',int2str(i));
end
Allrows=horzcat(P,S);
idx = find(cellfun(@isempty, Allrows));
Allrows(idx)={'none'};
Allrows=regexprep(Allrows,"'","");
AllColumnNames=horzcat(columnNamesP,columnNamesS);
SummaryTable=array2table(Allrows);
SummaryTable.Properties.VariableNames=AllColumnNames;
SummaryTable.Properties.RowNames=rownames;
save(strcat(saving_directory,modeldescription,'_summary.mat'),'SummaryTable')
writetable(SummaryTable,strcat(saving_directory,modeldescription,'_summary.csv'),'WriteRowNames',true);