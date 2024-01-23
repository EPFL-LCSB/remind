function elimIntermPMFile(saving_directory,corename)
% Eliminate intermediate files generated in PhenoMapping (recommended only
% after final result has been saved)
%
% USAGE:
%
%    elimIntermPMFile(saving_directory,corename)
%
% INPUT:
%    saving_directory:  Path to folder where results are saved
%    corename:          Name of the core module
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

if exist(strcat(saving_directory,corename,'_addGenes.mat'),'file') == 2
    delete(strcat(saving_directory,corename,'_addGenes.mat'))
end
if exist(strcat(saving_directory,corename,'_subsToGenes.mat'),'file') == 2
    delete(strcat(saving_directory,corename,'_subsToGenes.mat'))
end
if exist(strcat(saving_directory,corename,'_DPs.mat'),'file') == 2
    delete(strcat(saving_directory,corename,'_DPs.mat'))
end
if exist(strcat(saving_directory,corename,'_ess.mat'),'file') == 2
    delete(strcat(saving_directory,corename,'_ess.mat'))
end
if exist(strcat(saving_directory,'PhenoMappingDPMax_ess.mat'),'file') == 2
    delete(strcat(saving_directory,'PhenoMappingDPMax_ess.mat'))
end
if exist(strcat(saving_directory,'PhenoMappingDPMin_ess.mat'),'file') == 2
    delete(strcat(saving_directory,'PhenoMappingDPMin_ess.mat'))
end