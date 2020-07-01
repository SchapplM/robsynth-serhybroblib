% Calculate minimal parameter regressor of gravitation load for
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% taug_reg [1x9]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:15
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbarprisDE1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_gravloadJ_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisDE1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_gravloadJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:14:55
% EndTime: 2020-06-27 17:14:56
% DurationCPUTime: 0.11s
% Computational Cost: add. (257->38), mult. (210->52), div. (14->4), fcn. (28->2), ass. (0->31)
t19 = (qJ(1) ^ 2);
t21 = (pkin(3) ^ 2);
t45 = 2 * pkin(3) * qJ(1) + t19 + t21;
t14 = qJ(1) + pkin(3);
t35 = -pkin(2) + t14;
t36 = -pkin(2) - t14;
t26 = sqrt(-((pkin(1) + t36) * (pkin(1) + t35) * (pkin(1) - t35) * (pkin(1) - t36)));
t44 = g(1) * t26;
t43 = g(2) * t26;
t42 = 0.1e1 / pkin(1) / t26;
t22 = pkin(2) ^ 2;
t24 = pkin(1) ^ 2;
t39 = -t22 + t24;
t38 = t22 + t24;
t37 = 0.1e1 / (t14 ^ 2) * t42;
t34 = t14 / pkin(2) * t42;
t33 = t24 - t45;
t32 = -t37 / 0.2e1;
t7 = -t22 + t33;
t31 = t7 * t37 / 0.2e1;
t28 = t21 ^ 2;
t20 = pkin(3) * t21;
t18 = qJ(1) * t19;
t17 = pkin(1) - pkin(2);
t16 = pkin(1) + pkin(2);
t15 = -3 * t19;
t6 = t39 + t45;
t5 = t22 + t33;
t2 = g(1) * t6 + t43;
t1 = (g(2) * t6 - t44) * t31;
t3 = [0, t1, t2 * t7 * t32, t2 * t31, t1, (-(t20 + 4 * t21 * qJ(1) + (5 * t19 + t39) * pkin(3) + 2 * t18) * t44 + g(2) * (6 * t28 * qJ(1) + 2 * (7 * t19 - t38) * t20 + (-6 * t38 * qJ(1) + 16 * t18) * t21 - 2 * t18 * (-t19 + t38) + (t28 + (t17 ^ 2 + t15) * (t16 ^ 2 + t15)) * pkin(3))) * t32, 0, (g(2) * t5 + t44) * t34, -(g(1) * t5 - t43) * t34;];
taug_reg = t3;
