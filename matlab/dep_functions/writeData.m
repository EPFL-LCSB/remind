function writeData(filename,data,typeColumns,heading,typeHeading)
% Exports data to a file of the extension provided (.txt or .csv)
%
% USAGE:
%
%       writeData(filename,data,typeColumns,heading,typeHeading)
%
% INPUTS:
%    filename:        Filename with extension defined like .txt or .csv and
%                     path defined where you want to save it
%    data:            Data to export
%    typeColumns:     Indicate type of data in each column: string with %s,
%                     integer with %i, float with %f, ...
%
% OPTIONAL INPUTS:
%    heading:         Description of columns (default = empty / 
%                     no description)
%    typeHeading:     Indicate type of heading in each column: string with 
%                     %s, integer with %i, float with %f, ...
%
% OUTPUTS:
%    file:            Saved in path defined
%
% .. Author:
% Anush Chiappino 2018
% 

if (nargin < 4)
    heading = [];
end

if exist(filename,'file') == 2
    warning('the filename provided already exists. It will be overwritten')
end

fid = fopen(filename,'wt');
if fid>0
    if ~isempty(heading)
        fprintf(fid,strcat(typeHeading,'\n'),heading{1,:});
    end
    for k = 1:size(data,1)
        fprintf(fid,strcat(typeColumns,'\n'),data{k,:});
    end
    fclose(fid);
end
if exist(filename,'file') == 2
    fprintf('the data is now saved in a csv file \n');
else
    warning('the csv file was not created - check the saving_directory')
end