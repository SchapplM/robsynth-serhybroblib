% Convert vector of kinematic parameters to modified DH parameters of
% palh1m1TE
%
% Input:
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [16x1]
%   Rotation around z
% b_mdh [16x1]
%   Translation along z
% alpha_mdh [16x1]
%   Rotation around x
% a_mdh [16x1]
%   Translation along x
% theta_mdh [16x1]
%   Rotation around z
% d_mdh [16x1]
%   Translation along z
% qoffset_mdh [16x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = palh1m1TE_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(23,1)}
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_pkin2mdhparam: Kinematic parameters pkin have to be [23x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; -pkin(22); pkin(23); pkin(20); pkin(21); 0; 0; pi;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; -pkin(17); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; t1; 0; 0; t1; t1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(16); pkin(1); pkin(5); pkin(10); -pkin(15); pkin(1); 0; pkin(2); pkin(4); pkin(3); pkin(6); pkin(11); pkin(8); pkin(13); -pkin(9);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(14); 0; 0; 0; pkin(12); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; pi / 0.2e1; 0.3e1 / 0.2e1 * pi; 0; 0; 0; 0; 0; pi; pi; pi; pi; 0; 0; 0; 0;];
qoffset_mdh = t1;
