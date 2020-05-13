% Convert vector of kinematic parameters to modified DH parameters of
% palh2m2OL
%
% Input:
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [6x1]
%   Rotation around z
% b_mdh [6x1]
%   Translation along z
% alpha_mdh [6x1]
%   Rotation around x
% a_mdh [6x1]
%   Translation along x
% theta_mdh [6x1]
%   Rotation around z
% d_mdh [6x1]
%   Translation along z
% qoffset_mdh [6x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = palh2m2OL_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(5,1)}
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_pkin2mdhparam: Kinematic parameters pkin have to be [5x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = [0; pi / 0.2e1; 0; 0; 0; -pi / 0.2e1;];
alpha_mdh = t1;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(1); pkin(4); pkin(2); pkin(5); pkin(3);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
qoffset_mdh = t1;
