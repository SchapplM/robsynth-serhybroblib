% Starte die Modultests für alle generierten Robotermodelle
% 
% Moritz Schappler, schappler@imes.uni-hannover.de, 2018-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

clc
clear

%% Viergelenkkette
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar1', 'fourbar1_test.m'));

%% Fünfgelenkkette
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fivebar1', 'fivebar1_test.m'));

%% Palettierroboter mit drei geschlossenen Ketten (MPL800-Yaskawa)
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh1m1', 'palh1m1_test.m'));
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh1m2', 'palh1m2_test.m'));

%% Modellvergleich für Palettierroboter
test_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'test');
addpath(test_path);
run(fullfile(test_path, 'palletizer_mdlcomp'));

%% 2D-Delta
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'picker2Dm1', 'picker2Dm1_test.m'));

return
%% 3. Arm
% TODO: Hier gibt es noch Pfad-Probleme
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS5m5', 'testfcn'));
% Die Funktionen jedes Robotermodells werden für sich getestet
KAS5m5_test_everything_fixbase
KAS5m5u_test_everything_fixbase
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS5m5'));
KAS5m5_test_constr_dynamics

systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS5m6', 'testfcn'));
KAS5m6_test_everything_fixbase
KAS5m6u_test_everything_fixbase
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS5m6'));
KAS5m6_test_constr_dynamics

systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS5m7', 'testfcn'));
KAS7m1_test_everything_fixbase
KAS7m1u_test_everything_fixbase
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS5m7'));
KAS5m7_test_constr_dynamics

systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
cd(fullfile(systems_path, 'KAS7m1'));
KAS7m1_test_constr_dynamics

