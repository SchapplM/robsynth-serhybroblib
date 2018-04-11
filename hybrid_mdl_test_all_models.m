% Starte die Modultests für alle generierten Robotermodelle
% 
% Moritz Schappler, schappler@imes.uni-hannover.de, 2018-04
% (C) Institut für mechatronische Systeme, Universität Hannover

clc
clear

tb_path = fullfile(fileparts(which('hybrid_mdl_test_path_init.m')));
addpath(tb_path);
addpath(fullfile(tb_path, 'test'));
hybrid_mdl_test_path_init

%% Einzelne Modultests ausführen
% Die Funktionen jedes Robotermodells werden für sich getestet
KAS5m5_test_everything_fixbase
KAS5m5u_test_everything_fixbase
KAS5m6_test_everything_fixbase
KAS5m6u_test_everything_fixbase
KAS7m1_test_everything_fixbase
KAS7m1u_test_everything_fixbase

%% Tests zu den kinematischen Zwangsbedingungen ausführen
KAS5m5_test_constr_dynamics
KAS5m6_test_constr_dynamics
KAS5m7_test_constr_dynamics
KAS7m1_test_constr_dynamics