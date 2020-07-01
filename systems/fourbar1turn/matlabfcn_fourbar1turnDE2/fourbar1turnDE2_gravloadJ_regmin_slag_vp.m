% Calculate minimal parameter regressor of gravitation load for
% fourbar1turnDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% taug_reg [2x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:49
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnDE2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_gravloadJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:48:30
% EndTime: 2020-06-27 16:48:34
% DurationCPUTime: 0.50s
% Computational Cost: add. (4822->63), mult. (6790->127), div. (304->10), fcn. (1917->11), ass. (0->63)
t81 = pkin(4) ^ 2;
t45 = pkin(2) ^ 2;
t46 = pkin(1) ^ 2;
t38 = cos(qJ(2));
t68 = t38 * pkin(1);
t79 = -2 * pkin(2);
t62 = t68 * t79 + t46;
t31 = t45 + t62;
t29 = 0.1e1 / t31 ^ 2;
t80 = t29 / t81;
t28 = 0.1e1 / t31;
t61 = pkin(3) ^ 2 - t81;
t27 = t31 - t61;
t32 = -pkin(2) * t38 + pkin(1);
t36 = sin(qJ(2));
t76 = -pkin(3) - pkin(4);
t24 = (pkin(2) - t76) * (pkin(2) + t76) + t62;
t75 = -pkin(3) + pkin(4);
t25 = (pkin(2) - t75) * (pkin(2) + t75) + t62;
t48 = sqrt(-t24 * t25);
t65 = t36 * t48;
t16 = -pkin(2) * t65 + t27 * t32;
t70 = pkin(2) * t36;
t22 = t27 * t70;
t17 = t32 * t48 + t22;
t63 = t16 ^ 2 + t17 ^ 2;
t10 = t63 * t80;
t8 = t10 ^ (-0.1e1 / 0.2e1);
t78 = 0.2e1 * t36 ^ 2;
t64 = t48 * t38;
t69 = t36 * pkin(1);
t60 = pkin(2) * t69;
t67 = 0.1e1 / t48 * (-t24 - t25) * t60;
t3 = t22 + (-t64 + (0.2e1 * t32 * pkin(1) - t67) * t36) * pkin(2);
t4 = t32 * t67 + t45 * pkin(1) * t78 + (t27 * t38 + t65) * pkin(2);
t77 = 0.2e1 * (-0.2e1 * t28 * t60 * t63 + t16 * t3 + t17 * t4) * t8 / t10 * t80;
t39 = cos(qJ(1));
t74 = g(1) * t39;
t37 = sin(qJ(1));
t73 = g(2) * t37;
t72 = g(3) * t16;
t71 = g(3) * t17;
t44 = 0.1e1 / pkin(3);
t66 = t28 * t44;
t59 = t29 * t70;
t57 = pkin(1) * t8 * t59;
t56 = t73 + t74;
t55 = g(1) * t37 - g(2) * t39;
t54 = -t74 / 0.2e1 - t73 / 0.2e1;
t53 = 0.2e1 * t56;
t41 = 0.1e1 / pkin(4);
t52 = t55 * t8 * t41 * t28;
t33 = -pkin(2) + t68;
t26 = t31 + t61;
t23 = t26 * t69;
t18 = -t33 * t48 + t23;
t15 = -pkin(1) * t65 - t26 * t33;
t12 = 0.1e1 / t15 ^ 2;
t7 = qJ(2) + atan2(t18 * t66, t15 * t66);
t6 = cos(t7);
t5 = sin(t7);
t1 = 0.1e1 + (((t46 * pkin(2) * t78 - t33 * t67) * t28 + ((t26 * t38 + t65) * t28 - 0.2e1 * t18 * t59) * pkin(1)) / t15 - (t23 * t28 + (-t28 * t64 + ((t33 * t79 - t67) * t28 + t15 * t29 * t79) * t36) * pkin(1)) * t18 * t12) * pkin(3) / (t12 * t18 ^ 2 + 0.1e1) * t31 * t44;
t2 = [0, t55, t56, 0, 0, 0, 0, 0, t55 * t38, -t55 * t36, 0, 0, 0, 0, 0, 0, -t55 * t6, t55 * t5, 0, 0, 0, 0, 0, -t16 * t52, -t17 * t52; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t38 + t36 * t56, g(3) * t36 + t38 * t56, 0, 0, 0, 0, 0, 0, (g(3) * t6 - t5 * t56) * t1, (-g(3) * t5 - t56 * t6) * t1, 0, 0, 0, 0, 0, ((-t16 * t53 + 0.2e1 * t71) * t57 + ((-g(3) * t4 + t3 * t56) * t8 + (t71 / 0.2e1 + t54 * t16) * t77) * t28) * t41, ((-t17 * t53 - 0.2e1 * t72) * t57 + ((g(3) * t3 + t4 * t56) * t8 + (-t72 / 0.2e1 + t54 * t17) * t77) * t28) * t41;];
taug_reg = t2;
