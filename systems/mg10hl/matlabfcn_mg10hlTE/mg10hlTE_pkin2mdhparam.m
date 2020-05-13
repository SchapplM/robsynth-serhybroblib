% Convert vector of kinematic parameters to modified DH parameters of
% mg10hlTE
%
% Input:
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
%
% Output: MDH parameter in order of transformation
% beta_mdh [15x1]
%   Rotation around z
% b_mdh [15x1]
%   Translation along z
% alpha_mdh [15x1]
%   Rotation around x
% a_mdh [15x1]
%   Translation along x
% theta_mdh [15x1]
%   Rotation around z
% d_mdh [15x1]
%   Translation along z
% qoffset_mdh [15x1]
%   Offset on joint coordinate q

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:29
% Revision: f5729c120b58d1b0137ade1aa9321d1ea2b3cc0a (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = mg10hlTE_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(17,1)}
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlTE_pkin2mdhparam: Kinematic parameters pkin have to be [17x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t4 = -pi / 0.2e1;
t3 = pi / 0.2e1;
t1 = [0; 0; pkin(15); pkin(17); 0; t3; 0; t3; 0; 0; 0; 0; t4; 0; t4;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0; pkin(10); 0; 0; 0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t2 = -pi / 0.2e1;
t1 = pi / 0.2e1;
t3 = [0; t1; 0; 0; 0; t1; t2; t1; 0; t1; 0; 0; t2; 0; t2;];
alpha_mdh = t3;

% Aus parameters_mdh_a_matlab.m
t1 = [0; pkin(13); pkin(1); pkin(3); pkin(7); 0; 0; 0; 0; -pkin(14); pkin(6); pkin(5); 0; -pkin(4); 0;];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(11); 0; 0; 0; 0; pkin(12); 0; pkin(9); 0; 0; 0; 0; 0; 0; pkin(8);];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; pi / 0.2e1; pi; pi - pkin(16); pi; 0; -pi / 0.2e1; 0; 0; 0; pi; pi; 0; 0; 0;];
qoffset_mdh = t1;
