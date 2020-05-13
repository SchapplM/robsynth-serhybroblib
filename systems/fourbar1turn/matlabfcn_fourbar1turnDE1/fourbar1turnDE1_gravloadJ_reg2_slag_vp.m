% Calculate inertial parameters regressor of gravitation load for
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnDE1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t44 = cos(qJ(2));
t42 = sin(qJ(2));
t52 = pkin(2) ^ 2;
t53 = pkin(1) ^ 2;
t97 = pkin(2) * t44;
t84 = -0.2e1 * pkin(1) * t97 + t53;
t37 = t52 + t84;
t112 = pkin(3) ^ 2;
t113 = pkin(4) ^ 2;
t83 = t112 - t113;
t31 = t37 + t83;
t28 = pkin(1) * t42 * t31;
t39 = pkin(1) * t44 - pkin(2);
t107 = -pkin(3) - pkin(4);
t29 = (pkin(2) - t107) * (pkin(2) + t107) + t84;
t106 = pkin(4) - pkin(3);
t30 = (pkin(2) - t106) * (pkin(2) + t106) + t84;
t56 = sqrt(-t29 * t30);
t23 = -t39 * t56 + t28;
t87 = t44 * t56;
t98 = pkin(2) * t42;
t81 = pkin(1) * t98;
t90 = 0.1e1 / t56 * (-t29 - t30) * t81;
t6 = t28 + (-t87 + (-0.2e1 * t39 * pkin(2) - t90) * t42) * pkin(1);
t95 = t23 - t6;
t76 = t95 * t42;
t88 = t42 * t56;
t20 = -pkin(1) * t88 - t31 * t39;
t41 = t42 ^ 2;
t108 = 0.2e1 * t41;
t8 = -t39 * t90 + t53 * pkin(2) * t108 + (t31 * t44 + t88) * pkin(1);
t96 = t20 + t8;
t114 = t96 * t44 - t76;
t43 = sin(qJ(1));
t102 = g(2) * t43;
t45 = cos(qJ(1));
t103 = g(1) * t45;
t36 = t102 + t103;
t66 = 0.2e1 * t36;
t34 = 0.1e1 / t37 ^ 2;
t48 = 0.1e1 / t113;
t32 = t37 - t83;
t38 = pkin(1) - t97;
t21 = -pkin(2) * t88 + t32 * t38;
t27 = t32 * t98;
t22 = t38 * t56 + t27;
t85 = t21 ^ 2 + t22 ^ 2;
t14 = t85 * t48 * t34;
t10 = t14 ^ (-0.1e1 / 0.2e1);
t33 = 0.1e1 / t37;
t111 = t10 * t33;
t51 = 0.1e1 / t112;
t86 = t20 ^ 2 + t23 ^ 2;
t15 = t86 * t51 * t34;
t12 = t15 ^ (-0.1e1 / 0.2e1);
t110 = t12 * t33;
t75 = t34 * t81;
t50 = 0.1e1 / pkin(3);
t77 = t50 * t110;
t74 = t44 * t77;
t94 = t20 * t42;
t105 = (t23 * t74 + t77 * t94) * t45;
t104 = g(1) * t43;
t101 = g(2) * t45;
t100 = g(3) * t21;
t99 = g(3) * t22;
t93 = t20 * t44;
t92 = t23 * t42;
t91 = t23 * t44;
t89 = t42 * t44;
t82 = pkin(1) * pkin(2) * t34;
t80 = 0.2e1 * t34;
t7 = t27 + (-t87 + (0.2e1 * t38 * pkin(1) - t90) * t42) * pkin(2);
t72 = 0.4e1 * t33 * t75;
t9 = t38 * t90 + t52 * pkin(1) * t108 + (t32 * t44 + t88) * pkin(2);
t79 = ((t21 * t7 + t22 * t9) * t80 - t85 * t72) * t48 / t14 * t111;
t78 = 0.1e1 / t15 * ((t20 * t6 + t23 * t8) * t80 - t86 * t72) * t51 * t110;
t73 = -t101 + t104;
t71 = -t103 / 0.2e1 - t102 / 0.2e1;
t70 = t20 * t41 + t23 * t89;
t69 = t73 * t44;
t68 = t93 / 0.2e1 - t92 / 0.2e1;
t67 = -t91 / 0.2e1 - t94 / 0.2e1;
t65 = t96 * t42 + t95 * t44;
t47 = 0.1e1 / pkin(4);
t64 = t47 * t73 * t111;
t62 = 0.2e1 * t20 * t89 - 0.2e1 * t23 * t41;
t61 = -g(3) * t44 + t36 * t42;
t3 = t43 * t20 * t74;
t1 = [0, 0, 0, 0, 0, 0, t73, t36, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t73 * t42, -t36, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + (t92 * t104 - (t92 - t93) * t101) * t77, -g(2) * t105 - (-t91 - t94) * t77 * t104, -t36, pkin(2) * t69, 0, 0, 0, 0, 0, 0, -t21 * t64, -t22 * t64, -t36, t73 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, g(3) * t42 + t36 * t44, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t105 + ((g(3) * t67 - t36 * t68) * t78 + ((-0.2e1 * g(3) * t70 - t36 * t62) * t82 + (g(3) * t114 - (t42 * t8 - t44 * t6) * t103 - t65 * t102) * t33) * t12) * t50, -g(2) * t3 + ((-g(3) * t68 - t36 * t67) * t78 + ((-g(3) * t62 + t70 * t66) * t82 + (-g(3) * t65 - (t44 * t8 - t76) * t102 - t114 * t103) * t33) * t12) * t50, 0, t61 * pkin(2), 0, 0, 0, 0, 0, 0, ((t99 / 0.2e1 + t71 * t21) * t79 + ((-g(3) * t9 + t36 * t7) * t33 + (-t21 * t66 + 0.2e1 * t99) * t75) * t10) * t47, ((-t100 / 0.2e1 + t71 * t22) * t79 + ((g(3) * t7 + t36 * t9) * t33 + (-t22 * t66 - 0.2e1 * t100) * t75) * t10) * t47, 0, 0;];
taug_reg = t1;
