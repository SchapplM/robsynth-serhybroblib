% Pfade hinzufügen
% 
% For this file to work, it has to be executed as a complete file, not
% line-wise (because of the function `mfilename`)
% 
% Moritz Schappler, schappler@irt.uni-hannover.de, 2015-07
% (c) Institut für Regelungstechnik, Universität Hannover

%% Pfade hinzufügen

this_path = fileparts( mfilename('fullpath') );
addpath(this_path);

addpath(genpath(fullfile(this_path, 'tb_code')));
addpath(fullfile(this_path, 'test'));