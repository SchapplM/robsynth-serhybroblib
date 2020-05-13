% Jacobian of implicit kinematic constraints of
% fourbar2IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:37
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = fourbar2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:37:49
% EndTime: 2020-04-24 20:37:49
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->4), mult. (4->6), div. (5->2), fcn. (8->3), ass. (0->7)
t26 = -qJ(3) + qJ(1);
t28 = pkin(2) * sin(t26);
t22 = sin(qJ(2) + t26);
t21 = 0.1e1 / t22;
t27 = t21 / pkin(1);
t24 = sin(qJ(2));
t1 = [(-pkin(1) * t22 + t28) * t27; t24 * t21; (t24 * pkin(1) - t28) * t27;];
B21 = t1;
