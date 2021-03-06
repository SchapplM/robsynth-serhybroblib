% Convert vector of kinematic parameters to modified DH parameters of
% palh3m1OL
%
% Input:
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
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
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = palh3m1OL_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(16,1)}
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_pkin2mdhparam: Kinematic parameters pkin have to be [16x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; -pkin(15); 0; pi; pkin(16); pkin(14);];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; pkin(13); 0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; t1; 0; 0; t1; t1; 0; 0; 0; 0; 0; 0;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(12); pkin(1); pkin(4); pkin(8); -pkin(6); pkin(1); pkin(3); pkin(5); -pkin(7); pkin(2); pkin(9);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(11); 0; 0; 0; pkin(10); 0; 0; 0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; pi; 0; 0; 0; 0; pi; 0; 0; pi; 0;];
qoffset_mdh = t1;
