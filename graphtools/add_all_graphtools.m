function add_all_graphtools;
% main;
% 
% Initializes matlab paths to subfolders
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.

folder = pwd;

files = dir(folder);
for i=1:length(files)
    if files(i).isdir & strcmp(files(i).name,'.') == 0  && strcmp(files(i).name,'..') == 0
        addpath([folder '/' files(i).name]);
    end
end
