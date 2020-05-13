% Convert vector of modified DH parameters to kinematic parameter vector for
% mg10hlDE1
%
% Input:
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
%
% Output:
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 12:53
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = mg10hlDE1_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(3); NaN; a_mdh(4); -a_mdh(14); a_mdh(12); a_mdh(11); a_mdh(5); -d_mdh(15); d_mdh(8); b_mdh(7); d_mdh(1); d_mdh(6); a_mdh(2); -a_mdh(10); beta_mdh(3); -theta_mdh(4) - NaN + pi; beta_mdh(4);];
pkin = t1;
