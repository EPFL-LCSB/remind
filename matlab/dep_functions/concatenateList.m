function conCell = concatenateList(list, separator)
% Concatenate a list of elements in a cell array with the separator
% indicated
%
% USAGE:
%
%    conCell = concatenateList(list, separator)
%
% INPUT:
%    list:            List in a cell array
%    separator:       A character indicating the separator that we want to
%                     use to concatenate the list (e.g., '+' or ';' or ':' 
%                     or {' '},'+',{' '} etc.)
%
% OUTPUTS:
%    conCell:        concatenated list
%
% .. Author:
% Anush Chiappino-Pepe 2017
% 

conCell = {};
if (nargin < 2)
    separator = strcat(',',{' '});
end

for j = 1:length(list)
    if isempty(conCell)
        conCell = list{j};
    else
        conCell = strcat(conCell,separator,list{j});
    end
end