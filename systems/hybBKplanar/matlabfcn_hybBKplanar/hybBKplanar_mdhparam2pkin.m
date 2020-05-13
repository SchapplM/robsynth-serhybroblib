% Convert vector of modified DH parameters to kinematic parameter vector for
% hybBKplanar
%
% Input:
% beta_mdh [7x1]
%   Rotation around z
% b_mdh [7x1]
%   Translation along z
% alpha_mdh [7x1]
%   Rotation around x
% a_mdh [7x1]
%   Translation along x
% theta_mdh [7x1]
%   Rotation around z
% d_mdh [7x1]
%   Translation along z
% qoffset_mdh [7x1]
%   Offset on joint coordinate q
%
% Output:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,CF,ED]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 19:03
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = hybBKplanar_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(3); a_mdh(2); a_mdh(4); -a_mdh(6); a_mdh(5); a_mdh(7);];
pkin = t1;
