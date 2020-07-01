% Calculate minimal parameter regressor of gravitation load for
% fourbar1turnTE
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
% Datum: 2020-06-27 16:23
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnTE_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_gravloadJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:22:32
% EndTime: 2020-06-27 16:22:37
% DurationCPUTime: 0.42s
% Computational Cost: add. (1824->70), mult. (2632->141), div. (120->5), fcn. (791->6), ass. (0->68)
t30 = sin(qJ(1));
t65 = g(2) * t30;
t32 = cos(qJ(1));
t66 = g(1) * t32;
t46 = t65 + t66;
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t37 = pkin(2) ^ 2;
t38 = pkin(1) ^ 2;
t62 = pkin(2) * t31;
t55 = -0.2e1 * pkin(1) * t62 + t38;
t24 = t37 + t55;
t54 = pkin(3) ^ 2 - pkin(4) ^ 2;
t20 = t24 + t54;
t17 = pkin(1) * t29 * t20;
t26 = pkin(1) * t31 - pkin(2);
t71 = -pkin(3) - pkin(4);
t18 = (pkin(2) - t71) * (pkin(2) + t71) + t55;
t70 = -pkin(3) + pkin(4);
t19 = (pkin(2) - t70) * (pkin(2) + t70) + t55;
t39 = sqrt(-t18 * t19);
t56 = t31 * t39;
t63 = pkin(2) * t29;
t52 = pkin(1) * t63;
t60 = 0.1e1 / t39 * (-t18 - t19) * t52;
t1 = t17 + (-t56 + (-0.2e1 * t26 * pkin(2) - t60) * t29) * pkin(1);
t12 = -t26 * t39 + t17;
t74 = -t12 / 0.2e1;
t49 = t74 + t1 / 0.2e1;
t58 = t29 * t39;
t28 = t29 ^ 2;
t76 = 0.2e1 * t28;
t3 = -t26 * t60 + t38 * pkin(2) * t76 + (t31 * t20 + t58) * pkin(1);
t9 = -pkin(1) * t58 - t26 * t20;
t75 = t9 / 0.2e1;
t50 = t75 + t3 / 0.2e1;
t77 = t29 * t49 + t50 * t31;
t73 = t29 / 0.2e1;
t72 = t31 / 0.2e1;
t22 = 0.1e1 / t24;
t36 = 0.1e1 / pkin(3);
t59 = t22 * t36;
t47 = t59 * t72;
t51 = t29 * t59;
t69 = (t47 * t9 + t51 * t74) * t30;
t68 = (t12 * t47 + t51 * t75) * t32;
t67 = g(1) * t30;
t64 = g(2) * t32;
t61 = t31 * t9;
t57 = t31 * t12;
t23 = 0.1e1 / t24 ^ 2;
t53 = pkin(1) * pkin(2) * t23;
t48 = t23 * t52;
t45 = -t64 + t67;
t44 = -t12 * t28 + t29 * t61;
t43 = t28 * t9 + t29 * t57;
t42 = t66 / 0.2e1 + t65 / 0.2e1;
t34 = 0.1e1 / pkin(4);
t41 = (-t67 / 0.2e1 + t64 / 0.2e1) * t34 * t22;
t40 = t29 * t50 - t49 * t31;
t25 = pkin(1) - t62;
t21 = t24 - t54;
t16 = t21 * t63;
t11 = t25 * t39 + t16;
t10 = -pkin(2) * t58 + t25 * t21;
t4 = t25 * t60 + t37 * pkin(1) * t76 + (t31 * t21 + t58) * pkin(2);
t2 = t16 + (-t56 + (0.2e1 * t25 * pkin(1) - t60) * t29) * pkin(2);
t5 = [0, t45, t46, 0, 0, 0, 0, 0, t45 * t31, -t45 * t29, 0, 0, 0, 0, 0, 0, -g(1) * t69 - (-t61 / 0.2e1 + t12 * t73) * t59 * t64, -g(2) * t68 - (-t57 / 0.2e1 - t29 * t9 / 0.2e1) * t59 * t67, 0, 0, 0, 0, 0, t10 * t41, t11 * t41; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t31 + t29 * t46, g(3) * t29 + t46 * t31, 0, 0, 0, 0, 0, 0, -g(1) * t68 + ((-g(3) * t43 - t46 * t44) * t53 + (g(3) * t77 - (-t31 * t1 / 0.2e1 + t3 * t73) * t66 - t40 * t65) * t22) * t36, -g(2) * t69 + ((-g(3) * t44 + t46 * t43) * t53 + (-g(3) * t40 - (t1 * t73 + t3 * t72) * t65 - t77 * t66) * t22) * t36, 0, 0, 0, 0, 0, ((-g(3) * t4 / 0.2e1 + t42 * t2) * t22 + (g(3) * t11 - t10 * t46) * t48) * t34, ((g(3) * t2 / 0.2e1 + t42 * t4) * t22 + (-g(3) * t10 - t11 * t46) * t48) * t34;];
taug_reg = t5;
