% Calculate inertial parameters regressor of gravitation load for
% fourbar1DE1
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
% taug_reg [1x(1*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:57
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1DE1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1DE1_gravloadJ_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1DE1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1DE1_gravloadJ_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = cos(qJ(1));
t40 = pkin(2) * t25;
t37 = (-0.2e1 * t40 + pkin(1)) * pkin(1);
t45 = -pkin(3) - pkin(4);
t10 = (pkin(2) - t45) * (pkin(2) + t45) + t37;
t44 = -pkin(3) + pkin(4);
t11 = (pkin(2) - t44) * (pkin(2) + t44) + t37;
t33 = sqrt(-t10 * t11);
t24 = sin(qJ(1));
t41 = pkin(2) * t24;
t35 = pkin(1) * t41;
t47 = (-t10 - t11) * t35 / t33;
t21 = pkin(2) ^ 2 + t37;
t18 = 0.1e1 / t21;
t46 = t18 / 0.2e1;
t16 = -g(1) * t41 + g(2) * t40;
t12 = g(2) * pkin(1) - t16;
t34 = g(1) * t25 + g(2) * t24;
t17 = t34 * pkin(2);
t43 = -t12 * t47 - t17 * t33;
t13 = -g(1) * pkin(1) + t17;
t42 = -t13 * t47 - t16 * t33;
t39 = t12 * t18;
t38 = t13 * t18;
t36 = pkin(3) ^ 2 - pkin(4) ^ 2;
t30 = 0.1e1 / pkin(3);
t28 = 0.1e1 / pkin(4);
t20 = -g(1) * t24 + g(2) * t25;
t19 = 0.1e1 / t21 ^ 2;
t15 = t21 - t36;
t14 = t21 + t36;
t4 = t13 * t33;
t3 = t12 * t33;
t1 = [0, 0, 0, 0, 0, 0, -t20, t34, 0, 0, 0, 0, 0, 0, 0, 0, ((t14 * t16 + t43) * t46 + (t38 - (t14 * t13 - t3) * t19) * t35) * t30, ((-t14 * t17 + t42) * t46 + (-t39 - (-t14 * t12 - t4) * t19) * t35) * t30, 0, -pkin(2) * t20, 0, 0, 0, 0, 0, 0, ((-t15 * t16 + t43) * t46 + (-t38 - (-t15 * t13 - t3) * t19) * t35) * t28, ((t15 * t17 + t42) * t46 + (t39 - (t15 * t12 - t4) * t19) * t35) * t28, 0, 0;];
taug_reg = t1;
