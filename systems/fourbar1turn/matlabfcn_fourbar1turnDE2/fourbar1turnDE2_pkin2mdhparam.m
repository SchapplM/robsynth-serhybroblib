% Convert vector of kinematic parameters to modified DH parameters of
% fourbar1turnDE2
%
% Input:
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh] = fourbar1turnDE2_pkin2mdhparam(pkin)

%% Init
%#codegen
%$cgargs {zeros(5,1)}
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_pkin2mdhparam: Kinematic parameters pkin have to be [5x1] (double)');

%% Zuweisung der Parameter


% Aus parameters_mdh_beta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
beta_mdh = t1;

% Aus parameters_mdh_b_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
b_mdh = t1;

% Aus parameters_mdh_alpha_matlab.m
t1 = pi / 0.2e1;
t2 = [0; t1; 0; t1; 0; 0;];
alpha_mdh = t2;

% Aus parameters_mdh_a_matlab.m
t1 = [0; 0; pkin(2); pkin(1); pkin(3); pkin(4);];
a_mdh = t1;

% Aus parameters_mdh_theta_matlab.m
t1 = [0; 0; 0; 0; 0; 0;];
theta_mdh = t1;

% Aus parameters_mdh_d_matlab.m
t1 = [pkin(5); 0; 0; 0; 0; 0;];
d_mdh = t1;

% Aus parameters_mdh_qoffset_matlab.m
t1 = [0; 0; pi; 0; 0; 0;];
qoffset_mdh = t1;
