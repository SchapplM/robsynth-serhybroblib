% Convert vector of modified DH parameters to kinematic parameter vector for
% picker2Dm1OL
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = picker2Dm1OL_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2); a_mdh(3); a_mdh(4); a_mdh(10); a_mdh(5); a_mdh(9); a_mdh(7); beta_mdh(5);];
pkin = t1;
