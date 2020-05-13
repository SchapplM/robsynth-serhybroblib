% Convert vector of modified DH parameters to kinematic parameter vector for
% palh1m1OL
%
% Input:
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
%
% Output:
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = palh1m1OL_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(3); a_mdh(9); a_mdh(11); a_mdh(10); a_mdh(4); a_mdh(12); a_mdh(14); -a_mdh(16); a_mdh(5); a_mdh(13); d_mdh(5); a_mdh(15); d_mdh(1); -a_mdh(6); a_mdh(2); -b_mdh(6); beta_mdh(12); beta_mdh(13); -beta_mdh(10); beta_mdh(11);];
pkin = t1;
