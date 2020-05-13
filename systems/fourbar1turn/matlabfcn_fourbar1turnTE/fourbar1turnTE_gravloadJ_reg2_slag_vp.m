% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnTE_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnTE_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(qJ(1));
t67 = g(2) * t31;
t33 = cos(qJ(1));
t68 = g(1) * t33;
t24 = t67 + t68;
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t38 = pkin(2) ^ 2;
t39 = pkin(1) ^ 2;
t64 = pkin(2) * t32;
t57 = -0.2e1 * pkin(1) * t64 + t39;
t25 = t38 + t57;
t56 = pkin(3) ^ 2 - pkin(4) ^ 2;
t20 = t25 + t56;
t17 = pkin(1) * t30 * t20;
t27 = pkin(1) * t32 - pkin(2);
t73 = -pkin(3) - pkin(4);
t18 = (pkin(2) - t73) * (pkin(2) + t73) + t57;
t72 = -pkin(3) + pkin(4);
t19 = (pkin(2) - t72) * (pkin(2) + t72) + t57;
t40 = sqrt(-t18 * t19);
t58 = t32 * t40;
t65 = pkin(2) * t30;
t54 = pkin(1) * t65;
t62 = 0.1e1 / t40 * (-t18 - t19) * t54;
t1 = t17 + (-t58 + (-0.2e1 * t27 * pkin(2) - t62) * t30) * pkin(1);
t12 = -t27 * t40 + t17;
t76 = -t12 / 0.2e1;
t51 = t76 + t1 / 0.2e1;
t60 = t30 * t40;
t29 = t30 ^ 2;
t78 = 0.2e1 * t29;
t3 = -t27 * t62 + t39 * pkin(2) * t78 + (t32 * t20 + t60) * pkin(1);
t9 = -pkin(1) * t60 - t27 * t20;
t77 = t9 / 0.2e1;
t52 = t77 + t3 / 0.2e1;
t79 = t51 * t30 + t52 * t32;
t75 = t30 / 0.2e1;
t74 = t32 / 0.2e1;
t22 = 0.1e1 / t25;
t37 = 0.1e1 / pkin(3);
t61 = t22 * t37;
t49 = t61 * t74;
t53 = t30 * t61;
t71 = (t49 * t9 + t53 * t76) * t31;
t70 = (t12 * t49 + t53 * t77) * t33;
t69 = g(1) * t31;
t66 = g(2) * t33;
t63 = t32 * t9;
t59 = t32 * t12;
t23 = 0.1e1 / t25 ^ 2;
t55 = pkin(1) * pkin(2) * t23;
t50 = t23 * t54;
t48 = -t66 + t69;
t47 = -t12 * t29 + t30 * t63;
t46 = t29 * t9 + t30 * t59;
t45 = t68 / 0.2e1 + t67 / 0.2e1;
t44 = t48 * t32;
t35 = 0.1e1 / pkin(4);
t43 = (-t69 / 0.2e1 + t66 / 0.2e1) * t35 * t22;
t42 = -g(3) * t32 + t24 * t30;
t41 = t52 * t30 - t51 * t32;
t26 = pkin(1) - t64;
t21 = t25 - t56;
t16 = t21 * t65;
t11 = t26 * t40 + t16;
t10 = -pkin(2) * t60 + t26 * t21;
t4 = t26 * t62 + t38 * pkin(1) * t78 + (t32 * t21 + t60) * pkin(2);
t2 = t16 + (-t58 + (0.2e1 * t26 * pkin(1) - t62) * t30) * pkin(2);
t5 = [0, 0, 0, 0, 0, 0, t48, t24, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t48 * t30, -t24, 0, 0, 0, 0, 0, 0, 0, -g(1) * t71 - (-t63 / 0.2e1 + t12 * t75) * t61 * t66, -g(2) * t70 - (-t59 / 0.2e1 - t30 * t9 / 0.2e1) * t61 * t69, -t24, pkin(2) * t44, 0, 0, 0, 0, 0, 0, t10 * t43, t11 * t43, -t24, t48 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, g(3) * t30 + t24 * t32, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t70 + ((-g(3) * t46 - t24 * t47) * t55 + (g(3) * t79 - (-t32 * t1 / 0.2e1 + t3 * t75) * t68 - t41 * t67) * t22) * t37, -g(2) * t71 + ((-g(3) * t47 + t24 * t46) * t55 + (-g(3) * t41 - (t1 * t75 + t3 * t74) * t67 - t79 * t68) * t22) * t37, 0, t42 * pkin(2), 0, 0, 0, 0, 0, 0, ((-g(3) * t4 / 0.2e1 + t45 * t2) * t22 + (g(3) * t11 - t24 * t10) * t50) * t35, ((g(3) * t2 / 0.2e1 + t45 * t4) * t22 + (-g(3) * t10 - t24 * t11) * t50) * t35, 0, 0;];
taug_reg = t5;
