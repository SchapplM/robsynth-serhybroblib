% Convert vector of kinematic parameters to modified DH parameters of
% hybBKspatial
%
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED,L1,L2]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [10x1]
%   Rotation around z
% b_mdh [10x1]
%   Translation along z
% alpha_mdh [10x1]
%   Rotation around x
% a_mdh [10x1]
%   Translation along x
% theta_mdh [10x1]
%   Rotation around z
% d_mdh [10x1]
%   Translation along z
% qoffset_mdh [10x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 19:31
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = hybBKspatial_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(7,1)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'hybBKspatial_pkin2mdhparam: Kinematic parameters pkin have to be [7x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; -pi / 0.2e1; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; t1; 0; t1; 0; 0; -pi / 0.2e1; t1; 0; 0;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; -pkin(1) / 0.2e1; pkin(2); pkin(1) / 0.2e1; pkin(3); pkin(7); 0; 0; -pkin(4); pkin(5);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(6); 0; 0; 0; 0; 0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; pi; 0; pi; 0; 0; 0; 0; 0;];
qoffset_mdh = t1;
