% Convert vector of modified DH parameters to kinematic parameter vector for
% picker2Dm2DE2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05,phi1]';

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:02
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function pkin = picker2Dm2DE2_mdhparam2pkin(beta_mdh, b_mdh, alpha_mdh, a_mdh, theta_mdh, d_mdh, qoffset_mdh)

% Aus parameter_kin_from_mdh_matlab.m
t1 = [a_mdh(2); a_mdh(3); a_mdh(4); a_mdh(10); a_mdh(5); a_mdh(9); a_mdh(7); beta_mdh(5); NaN;];
pkin = t1;
