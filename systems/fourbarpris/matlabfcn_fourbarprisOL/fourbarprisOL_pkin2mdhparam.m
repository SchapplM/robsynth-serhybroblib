% Convert vector of kinematic parameters to modified DH parameters of
% fourbarprisOL
%
% Input:
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [5x1]
%   Rotation around z
% b_mdh [5x1]
%   Translation along z
% alpha_mdh [5x1]
%   Rotation around x
% a_mdh [5x1]
%   Translation along x
% theta_mdh [5x1]
%   Rotation around z
% d_mdh [5x1]
%   Translation along z
% qoffset_mdh [5x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = fourbarprisOL_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(3,1)}
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_pkin2mdhparam: Kinematic parameters pkin have to be [3x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; -pi / 0.2e1; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; pkin(3);];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = [0; -pi / 0.2e1; 0; 0; pi / 0.2e1;];
alpha_mdh = t1;

% Aus parameters_mdh_a_matlab.m
t1 = [pkin(1); 0; 0; pkin(2); 0;];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; pi;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [pi; 0; pi; pi; pi / 0.2e1;];
qoffset_mdh = t1;
