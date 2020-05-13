% Calculate inertial parameters regressor of gravitation load for
% palh1m2TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = palh1m2TE_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2TE_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_gravloadJ_reg2_slag_vp: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t78 = sin(qJ(1));
t84 = cos(qJ(1));
t143 = -g(1) * t78 + g(2) * t84;
t117 = sin(pkin(20));
t71 = sin(pkin(21));
t73 = cos(pkin(21));
t74 = cos(pkin(20));
t46 = -t117 * t73 - t74 * t71;
t47 = -t117 * t71 + t74 * t73;
t70 = sin(pkin(22));
t72 = cos(pkin(22));
t112 = t46 * t70 + t47 * t72;
t113 = -t46 * t72 + t70 * t47;
t79 = sin(pkin(18));
t85 = cos(pkin(18));
t146 = -t112 * t79 + t113 * t85;
t153 = t143 * t146;
t57 = g(1) * t84 + g(2) * t78;
t80 = sin(pkin(17));
t86 = cos(pkin(17));
t50 = t79 * t86 - t85 * t80;
t51 = t80 * t79 + t86 * t85;
t77 = sin(qJ(2));
t83 = cos(qJ(2));
t150 = t50 * t77 + t51 * t83;
t118 = sin(pkin(19));
t119 = cos(pkin(19));
t76 = sin(qJ(3));
t82 = cos(qJ(3));
t48 = t82 * t118 + t76 * t119;
t49 = -t76 * t118 + t82 * t119;
t149 = -t48 * t77 + t49 * t83;
t148 = t48 * t83 + t49 * t77;
t18 = t148 * g(1);
t19 = t148 * g(2);
t142 = -g(3) * t149 + t18 * t84 + t19 * t78;
t141 = g(3) * t146;
t140 = pkin(1) * t77;
t66 = pkin(5) * t76 + pkin(1);
t124 = t66 * t77;
t123 = t66 * t83;
t121 = t77 * t82;
t120 = t82 * t83;
t116 = pkin(5) * t121;
t115 = pkin(5) * t120;
t101 = t112 * t85 + t113 * t79;
t1 = t101 * g(1);
t75 = sin(qJ(4));
t81 = cos(qJ(4));
t114 = g(2) * t81 - t1 * t75;
t109 = -pkin(15) - t115;
t108 = t150 * g(1);
t107 = t150 * g(2);
t2 = t101 * g(2);
t105 = -g(1) * t81 - t2 * t75;
t104 = g(1) * t75 - t2 * t81;
t103 = g(2) * t75 + t1 * t81;
t38 = t48 * pkin(2);
t39 = t49 * pkin(2);
t99 = t38 * t77 - t39 * t83;
t97 = t50 * t83 - t51 * t77;
t94 = t76 * t83 + t121;
t93 = t76 * t77 - t120;
t27 = g(3) * t77 + t57 * t83;
t88 = t99 * g(3) + t57 * (t38 * t83 + t39 * t77);
t87 = t148 * g(3) + t57 * t149;
t68 = g(1) * t140;
t67 = g(2) * t140;
t59 = g(1) * t116;
t58 = g(2) * t116;
t56 = pkin(9) * t74 - pkin(11) * t117;
t55 = t117 * pkin(9) + pkin(11) * t74;
t54 = g(2) * t124;
t53 = g(1) * t124;
t52 = t143 * pkin(15);
t45 = t143 * t83;
t44 = t143 * t77;
t43 = t94 * g(2);
t42 = t93 * g(2);
t41 = t94 * g(1);
t40 = t93 * g(1);
t26 = g(3) * t83 - t57 * t77;
t23 = t27 * pkin(1);
t22 = -(-pkin(15) * g(1) + t68) * t78 + (-g(2) * pkin(15) + t67) * t84;
t21 = -t55 * t71 + t56 * t73;
t20 = t55 * t73 + t56 * t71;
t11 = -t41 * t78 + t43 * t84;
t10 = -t40 * t78 + t42 * t84;
t9 = t94 * g(3) - t40 * t84 - t42 * t78;
t8 = t93 * g(3) + t41 * t84 + t43 * t78;
t7 = t58 * t78 + t59 * t84 + (-g(3) * t120 + t27 * t76) * pkin(5);
t6 = (g(1) * t123 + t59) * t84 + (g(2) * t123 + t58) * t78 + (-t115 + t124) * g(3);
t3 = [0, 0, 0, 0, 0, 0, -t143, t57, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, -t57, -t52, 0, 0, 0, 0, 0, 0, t10, t11, -t57, t22, 0, 0, 0, 0, 0, 0, -t1 * t78 + t2 * t84, t153, -t57, -(t109 * g(1) + t53) * t78 + (t109 * g(2) + t54) * t84, 0, 0, 0, 0, 0, 0, -t103 * t78 - t104 * t84, t105 * t84 - t114 * t78, -t153, -t53 * t78 + t54 * t84 + t143 * ((t20 * t72 + t21 * t70) * t79 + (-t20 * t70 + t21 * t72) * t85 + t109), 0, 0, 0, 0, 0, 0, -t143 * t97, t107 * t84 - t108 * t78, -t57, t143 * pkin(14), 0, 0, 0, 0, 0, 0, t143 * (t70 * t79 + t72 * t85), t143 * (t70 * t85 - t72 * t79), -t57, t22, 0, 0, 0, 0, 0, 0, -t143 * t149, -t18 * t78 + t19 * t84, -t57, -t52, 0, 0, 0, 0, 0, 0, t44, t45, -t57, t143 * (-pkin(15) + t99), 0, 0, 0, 0, 0, 0, t10, t11, -t57, t67 * t84 - t68 * t78 + t143 * (-pkin(15) + ((t70 * t73 + t71 * t72) * t79 + (-t70 * t71 + t72 * t73) * t85) * pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, -t97 * g(3) + t107 * t78 + t108 * t84, t150 * g(3) + t57 * t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, t142, t87, 0, 0, 0, 0, 0, 0, 0, 0, t27, t26, 0, t88, 0, 0, 0, 0, 0, 0, t8, t9, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t78 + t114 * t84 + t75 * t141, -t103 * t84 + t104 * t78 + t81 * t141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
taug_reg = t3;
