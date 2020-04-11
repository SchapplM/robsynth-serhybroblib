% Starte die Modultests für alle generierten Robotermodelle
% 
% Moritz Schappler, schappler@imes.uni-hannover.de, 2018-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

%% Viergelenkkette
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar1', 'fourbar1_test.m'));
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar2', 'fourbar2_test.m'));
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbarpris', 'fourbarpris_test.m'));
%% Fünfgelenkkette
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fivebar1', 'fivebar1_test.m'));

%% Minimalbeispiele und Testsysteme
% Als Teil der größeren Roboter zum Testen und zur Entwicklung
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'fourbar1turn', 'fourbar1turn_test.m'));

%% Palettierroboter mit drei geschlossenen Ketten (Typ Yaskawa MPL800)
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh1m1', 'palh1m1_test.m'));
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh1m2', 'palh1m2_test.m'));

%% Palettierroboter mit zwei geschlossenen Ketten (Typ KUKA KR-700PA)
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh3m1', 'palh3m1_test.m'));
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'palh3m2', 'palh3m2_test.m'));

%% Modellvergleich für Palettierroboter
test_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'test');
addpath(test_path);
run(fullfile(test_path, 'palletizer_mdlcomp'));

%% 2D-Palettierroboter
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'picker2Dm1', 'picker2Dm1_test.m'));
systems_path = fullfile(fileparts(which('hybroblib_path_init.m')), 'systems');
run(fullfile(systems_path, 'picker2Dm2', 'picker2Dm2_test.m'));
