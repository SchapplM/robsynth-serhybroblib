tb_path = fullfile(fileparts(which('hybrid_mdl_test_path_init.m')));
addpath(tb_path);
hybrid_mdl_test_path_init

KAS_path = fullfile(fileparts(which('KAS_path.m')));
rmpath(genpath(KAS_path));

KAS5m5_test_everything_fixbase
KAS5m5u_test_everything_fixbase
KAS5m6_test_everything_fixbase
KAS5m6u_test_everything_fixbase
KAS7m1_test_everything_fixbase
KAS7m1u_test_everything_fixbase

