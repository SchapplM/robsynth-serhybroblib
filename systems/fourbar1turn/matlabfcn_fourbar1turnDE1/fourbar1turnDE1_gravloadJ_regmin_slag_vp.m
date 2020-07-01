% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [2x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 16:36
% Revision: 75f93b5b4b0ac6379b75b4546e5e7b5b01e11d8f (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = fourbar1turnDE1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_gravloadJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1turnDE1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 16:35:23
% EndTime: 2020-06-27 16:35:29
% DurationCPUTime: 0.90s
% Computational Cost: add. (8832->87), mult. (12850->174), div. (612->9), fcn. (3605->10), ass. (0->88)
t43 = cos(qJ(2));
t41 = sin(qJ(2));
t51 = pkin(2) ^ 2;
t52 = pkin(1) ^ 2;
t95 = pkin(2) * t43;
t82 = -0.2e1 * pkin(1) * t95 + t52;
t36 = t51 + t82;
t110 = pkin(3) ^ 2;
t111 = pkin(4) ^ 2;
t81 = t110 - t111;
t31 = t36 + t81;
t28 = pkin(1) * t41 * t31;
t38 = pkin(1) * t43 - pkin(2);
t105 = -pkin(3) - pkin(4);
t29 = (pkin(2) - t105) * (pkin(2) + t105) + t82;
t104 = pkin(4) - pkin(3);
t30 = (pkin(2) - t104) * (pkin(2) + t104) + t82;
t55 = sqrt(-t29 * t30);
t23 = -t38 * t55 + t28;
t85 = t43 * t55;
t96 = pkin(2) * t41;
t79 = pkin(1) * t96;
t88 = 0.1e1 / t55 * (-t29 - t30) * t79;
t6 = t28 + (-t85 + (-0.2e1 * t38 * pkin(2) - t88) * t41) * pkin(1);
t93 = t23 - t6;
t74 = t93 * t41;
t86 = t41 * t55;
t20 = -pkin(1) * t86 - t31 * t38;
t40 = t41 ^ 2;
t106 = 0.2e1 * t40;
t8 = -t38 * t88 + t52 * pkin(2) * t106 + (t31 * t43 + t86) * pkin(1);
t94 = t20 + t8;
t112 = t94 * t43 - t74;
t42 = sin(qJ(1));
t100 = g(2) * t42;
t44 = cos(qJ(1));
t101 = g(1) * t44;
t71 = t100 + t101;
t64 = 0.2e1 * t71;
t34 = 0.1e1 / t36 ^ 2;
t47 = 0.1e1 / t111;
t32 = t36 - t81;
t37 = pkin(1) - t95;
t21 = -pkin(2) * t86 + t32 * t37;
t27 = t32 * t96;
t22 = t37 * t55 + t27;
t83 = t21 ^ 2 + t22 ^ 2;
t14 = t83 * t47 * t34;
t10 = t14 ^ (-0.1e1 / 0.2e1);
t33 = 0.1e1 / t36;
t109 = t10 * t33;
t50 = 0.1e1 / t110;
t84 = t20 ^ 2 + t23 ^ 2;
t15 = t84 * t50 * t34;
t12 = t15 ^ (-0.1e1 / 0.2e1);
t108 = t12 * t33;
t73 = t34 * t79;
t49 = 0.1e1 / pkin(3);
t75 = t49 * t108;
t72 = t43 * t75;
t92 = t20 * t41;
t103 = (t23 * t72 + t75 * t92) * t44;
t102 = g(1) * t42;
t99 = g(2) * t44;
t98 = g(3) * t21;
t97 = g(3) * t22;
t91 = t20 * t43;
t90 = t23 * t41;
t89 = t23 * t43;
t87 = t41 * t43;
t80 = pkin(1) * pkin(2) * t34;
t78 = 0.2e1 * t34;
t69 = 0.4e1 * t33 * t73;
t7 = t27 + (-t85 + (0.2e1 * t37 * pkin(1) - t88) * t41) * pkin(2);
t9 = t37 * t88 + t51 * pkin(1) * t106 + (t32 * t43 + t86) * pkin(2);
t77 = ((t21 * t7 + t22 * t9) * t78 - t83 * t69) * t47 / t14 * t109;
t76 = 0.1e1 / t15 * ((t20 * t6 + t23 * t8) * t78 - t84 * t69) * t50 * t108;
t70 = -t99 + t102;
t68 = -t101 / 0.2e1 - t100 / 0.2e1;
t67 = t20 * t40 + t23 * t87;
t66 = t91 / 0.2e1 - t90 / 0.2e1;
t65 = -t89 / 0.2e1 - t92 / 0.2e1;
t63 = t94 * t41 + t93 * t43;
t46 = 0.1e1 / pkin(4);
t62 = t46 * t70 * t109;
t60 = 0.2e1 * t20 * t87 - 0.2e1 * t23 * t40;
t3 = t42 * t20 * t72;
t1 = [0, t70, t71, 0, 0, 0, 0, 0, t70 * t43, -t70 * t41, 0, 0, 0, 0, 0, 0, -g(1) * t3 + (t90 * t102 - (t90 - t91) * t99) * t75, -g(2) * t103 - (-t89 - t92) * t75 * t102, 0, 0, 0, 0, 0, -t21 * t62, -t22 * t62; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t43 + t71 * t41, g(3) * t41 + t71 * t43, 0, 0, 0, 0, 0, 0, -g(1) * t103 + ((g(3) * t65 - t71 * t66) * t76 + ((-0.2e1 * g(3) * t67 - t71 * t60) * t80 + (g(3) * t112 - (t41 * t8 - t43 * t6) * t101 - t63 * t100) * t33) * t12) * t49, -g(2) * t3 + ((-g(3) * t66 - t71 * t65) * t76 + ((-g(3) * t60 + t67 * t64) * t80 + (-g(3) * t63 - (t43 * t8 - t74) * t100 - t112 * t101) * t33) * t12) * t49, 0, 0, 0, 0, 0, ((t97 / 0.2e1 + t68 * t21) * t77 + ((-g(3) * t9 + t71 * t7) * t33 + (-t21 * t64 + 0.2e1 * t97) * t73) * t10) * t46, ((-t98 / 0.2e1 + t68 * t22) * t77 + ((g(3) * t7 + t71 * t9) * t33 + (-t22 * t64 - 0.2e1 * t98) * t73) * t10) * t46;];
taug_reg = t1;
