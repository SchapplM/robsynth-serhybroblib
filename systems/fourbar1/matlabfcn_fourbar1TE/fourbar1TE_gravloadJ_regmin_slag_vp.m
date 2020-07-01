% Calculate minimal parameter regressor of gravitation load for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% taug_reg [1x9]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1TE_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_gravloadJ_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1TE_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:21:05
% EndTime: 2020-06-26 17:21:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (318->38), mult. (492->56), div. (16->5), fcn. (140->4), ass. (0->33)
t24 = cos(qJ(1));
t39 = pkin(2) * t24;
t36 = (-0.2e1 * t39 + pkin(1)) * pkin(1);
t44 = -pkin(3) - pkin(4);
t10 = (pkin(2) - t44) * (pkin(2) + t44) + t36;
t43 = -pkin(3) + pkin(4);
t11 = (pkin(2) - t43) * (pkin(2) + t43) + t36;
t32 = sqrt(-t10 * t11);
t23 = sin(qJ(1));
t40 = pkin(2) * t23;
t34 = pkin(1) * t40;
t46 = (-t10 - t11) * t34 / t32;
t20 = pkin(2) ^ 2 + t36;
t18 = 0.1e1 / t20;
t45 = t18 / 0.2e1;
t16 = -g(1) * t40 + g(2) * t39;
t12 = g(2) * pkin(1) - t16;
t33 = g(1) * t24 + g(2) * t23;
t17 = t33 * pkin(2);
t42 = -t12 * t46 - t17 * t32;
t13 = -g(1) * pkin(1) + t17;
t41 = -t13 * t46 - t16 * t32;
t38 = t12 * t18;
t37 = t13 * t18;
t35 = pkin(3) ^ 2 - pkin(4) ^ 2;
t29 = 0.1e1 / pkin(3);
t27 = 0.1e1 / pkin(4);
t19 = 0.1e1 / t20 ^ 2;
t15 = t20 - t35;
t14 = t20 + t35;
t4 = t13 * t32;
t3 = t12 * t32;
t1 = [0, g(1) * t23 - g(2) * t24, t33, 0, ((t14 * t16 + t42) * t45 + (t37 - (t14 * t13 - t3) * t19) * t34) * t29, ((-t14 * t17 + t41) * t45 + (-t38 - (-t14 * t12 - t4) * t19) * t34) * t29, 0, ((-t15 * t16 + t42) * t45 + (-t37 - (-t15 * t13 - t3) * t19) * t34) * t27, ((t15 * t17 + t41) * t45 + (t38 - (t15 * t12 - t4) * t19) * t34) * t27;];
taug_reg = t1;
