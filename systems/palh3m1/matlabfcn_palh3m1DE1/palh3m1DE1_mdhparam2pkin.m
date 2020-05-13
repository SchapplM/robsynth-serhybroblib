% Convert vector of modified DH parameters to kinematic parameter vector for
% palh3m1DE1
%
% Input:
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
%
% Output:
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = palh3m1DE1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(3); a_mdh(11); a_mdh(8); a_mdh(4); NaN; a_mdh(9); -a_mdh(6); -a_mdh(10); a_mdh(5); a_mdh(12); d_mdh(5); d_mdh(1); a_mdh(2); b_mdh(6); NaN; NaN; beta_mdh(12); -beta_mdh(8); beta_mdh(11);];
pkin = t1;
