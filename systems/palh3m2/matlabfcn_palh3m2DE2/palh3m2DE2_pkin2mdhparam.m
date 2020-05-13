% Convert vector of kinematic parameters to modified DH parameters of
% palh3m2DE2
%
% Input:
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [12x1]
%   Rotation around z
% b_mdh [12x1]
%   Translation along z
% alpha_mdh [12x1]
%   Rotation around x
% a_mdh [12x1]
%   Translation along x
% theta_mdh [12x1]
%   Rotation around z
% d_mdh [12x1]
%   Translation along z
% qoffset_mdh [12x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = palh3m2DE2_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(18,1)}
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_pkin2mdhparam: Kinematic parameters pkin have to be [18x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; -pkin(17); pkin(18); pkin(16); 0; pi;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; pkin(13); 0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; t1; 0; 0; t1; t1; 0; 0; 0; 0; 0; 0;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(12); pkin(1); pkin(4); pkin(8); -pkin(6); pkin(1); pkin(3); pkin(2); pkin(9); pkin(5); -pkin(7);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(11); 0; 0; 0; pkin(10); 0; 0; 0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; pi; 0; 0; 0; 0; pi; pi; 0; 0; 0;];
qoffset_mdh = t1;
