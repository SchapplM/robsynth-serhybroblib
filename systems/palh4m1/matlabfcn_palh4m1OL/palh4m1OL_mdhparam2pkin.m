% Convert vector of modified DH parameters to kinematic parameter vector for
% palh4m1OL
%
% Input:
% beta_mdh [9x1]
%   Rotation around z
% b_mdh [9x1]
%   Translation along z
% alpha_mdh [9x1]
%   Rotation around x
% a_mdh [9x1]
%   Translation along x
% theta_mdh [9x1]
%   Rotation around z
% d_mdh [9x1]
%   Translation along z
% qoffset_mdh [9x1]
%   Offset on joint coordinate q
%
% Output:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:04
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = palh4m1OL_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(9); a_mdh(8); a_mdh(5); d_mdh(6); d_mdh(1); a_mdh(7); -a_mdh(2);];
pkin = t1;
