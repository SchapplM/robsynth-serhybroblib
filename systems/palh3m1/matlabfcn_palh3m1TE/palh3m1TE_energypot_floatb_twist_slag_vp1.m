% Calculate potential energy for
% palh3m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [19x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DA,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1TE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:18:49
% EndTime: 2020-04-17 15:18:57
% DurationCPUTime: 4.97s
% Computational Cost: add. (78772->185), mult. (118428->266), div. (5616->6), fcn. (74977->24), ass. (0->107)
t124 = -pkin(11) - rSges(6,3);
t131 = 0.1e1 / pkin(2);
t130 = pkin(2) - pkin(6);
t129 = -pkin(6) - pkin(2);
t128 = -pkin(8) - pkin(10);
t127 = -pkin(8) + pkin(10);
t59 = cos(qJ(2));
t126 = pkin(1) * t59;
t54 = sin(qJ(2));
t56 = sin(pkin(16));
t61 = cos(pkin(16));
t42 = t54 * t56 - t59 * t61;
t125 = pkin(5) * t42;
t123 = sin(pkin(15));
t120 = (-0.2e1 * t125 + pkin(1)) * pkin(1);
t109 = pkin(5) ^ 2 + t120;
t31 = 0.1e1 / t109;
t122 = t31 / pkin(6);
t29 = sqrt(-((pkin(5) - t130) * (pkin(5) + t130) + t120) * ((pkin(5) - t129) * (pkin(5) + t129) + t120));
t44 = t54 * t61 + t56 * t59;
t121 = t44 * t29;
t119 = -pkin(2) ^ 2 + pkin(6) ^ 2;
t118 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t117 = cos(pkin(19));
t116 = sin(pkin(19));
t115 = cos(pkin(17));
t114 = sin(pkin(17));
t113 = pkin(18) + pkin(19);
t112 = pkin(12) + r_base(3);
t55 = sin(qJ(1));
t111 = t55 * pkin(13) + r_base(2);
t60 = cos(qJ(1));
t110 = t60 * pkin(13) + r_base(1);
t108 = cos(pkin(15)) / 0.2e1;
t107 = pkin(1) - t125;
t106 = t54 * pkin(1) + t112;
t105 = t55 * t126 + t111;
t104 = t60 * t126 + t110;
t103 = cos(t113);
t102 = sin(t113);
t53 = sin(qJ(3));
t58 = cos(qJ(3));
t96 = t53 * t59 + t54 * t58;
t101 = -pkin(4) * t96 + t106;
t41 = t53 * t54 - t58 * t59;
t35 = t41 * t55;
t100 = t35 * pkin(4) + t105;
t37 = t41 * t60;
t99 = t37 * pkin(4) + t104;
t98 = t109 - t119;
t97 = rSges(3,1) * t59 - rSges(3,2) * t54;
t92 = t131 * (-pkin(5) * t121 + t107 * t98);
t90 = -t92 / 0.2e1;
t91 = t131 * (pkin(5) * t44 * t98 + t107 * t29) / 0.2e1;
t21 = (t116 * t91 + t117 * t90) * t31;
t89 = t92 / 0.2e1;
t22 = (t116 * t89 + t117 * t91) * t31;
t20 = t21 * t59 - t22 * t54;
t19 = t21 * t54 + t22 * t59;
t30 = t109 + t119;
t38 = pkin(1) * t42 - pkin(5);
t95 = pkin(1) * t30 * t44 - t29 * t38;
t94 = -pkin(1) * t121 - t30 * t38;
t23 = (t94 * t108 + t95 * t123 / 0.2e1) * t122;
t24 = (t95 * t108 - t94 * t123 / 0.2e1) * t122;
t93 = rSges(7,1) * t23 - rSges(7,2) * t24 - pkin(7);
t88 = t53 * t89 + t58 * t91;
t87 = t53 * t91 + t58 * t90;
t86 = t31 * (-t102 * t88 - t103 * t87);
t85 = t31 * (t102 * t87 - t103 * t88);
t84 = pkin(4) * t86;
t83 = pkin(3) - t84;
t82 = -pkin(3) * t86 + pkin(4);
t81 = -0.2e1 * pkin(3) * t84 + pkin(4) ^ 2;
t80 = pkin(3) ^ 2 + t81;
t79 = 0.1e1 / t80;
t78 = 0.1e1 / pkin(8) * t79;
t77 = 0.1e1 / pkin(10) * t79;
t76 = t80 + t118;
t75 = t80 - t118;
t74 = sqrt(-((pkin(3) - t127) * (pkin(3) + t127) + t81) * ((pkin(3) - t128) * (pkin(3) + t128) + t81));
t73 = t74 * t85;
t72 = (-pkin(4) * t73 + t83 * t75) * t78;
t71 = (-pkin(3) * t73 + t82 * t76) * t77;
t70 = -(pkin(4) * t75 * t85 + t83 * t74) * t78 / 0.2e1;
t69 = (pkin(3) * t76 * t85 + t82 * t74) * t77 / 0.2e1;
t57 = cos(qJ(4));
t52 = sin(qJ(4));
t51 = cos(pkin(18));
t50 = sin(pkin(18));
t36 = t96 * t60;
t34 = t96 * t55;
t18 = t19 * t60;
t17 = t20 * t60;
t16 = t19 * t55;
t15 = t20 * t55;
t10 = -t51 * t72 / 0.2e1 + t50 * t70;
t9 = t50 * t72 / 0.2e1 + t51 * t70;
t8 = t115 * t69 + t114 * t71 / 0.2e1;
t7 = -t115 * t71 / 0.2e1 + t114 * t69;
t6 = t41 * t7 + t8 * t96;
t5 = t41 * t8 - t7 * t96;
t4 = t36 * t7 - t37 * t8;
t3 = t36 * t8 + t37 * t7;
t2 = t34 * t7 - t35 * t8;
t1 = t34 * t8 + t35 * t7;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t60 - rSges(2,2) * t55 + r_base(1)) + g(2) * (rSges(2,1) * t55 + rSges(2,2) * t60 + r_base(2)) + g(3) * (rSges(2,3) + t112)) - m(3) * (g(1) * (rSges(3,3) * t55 + t97 * t60 + t110) + g(2) * (-rSges(3,3) * t60 + t97 * t55 + t111) + g(3) * (rSges(3,1) * t54 + rSges(3,2) * t59 + t112)) - m(4) * (g(1) * (rSges(4,1) * t37 + rSges(4,2) * t36 + rSges(4,3) * t55 + t104) + g(2) * (rSges(4,1) * t35 + rSges(4,2) * t34 - rSges(4,3) * t60 + t105) + g(3) * (-rSges(4,1) * t96 + rSges(4,2) * t41 + t106)) - m(5) * (g(1) * (rSges(5,1) * t3 + rSges(5,2) * t4 + rSges(5,3) * t55 + t99) + g(2) * (rSges(5,1) * t1 + rSges(5,2) * t2 - rSges(5,3) * t60 + t100) + g(3) * (rSges(5,1) * t5 + rSges(5,2) * t6 + t101)) - m(6) * (g(1) * (t3 * pkin(9) + (t3 * t57 + t52 * t55) * rSges(6,1) + (-t3 * t52 + t55 * t57) * rSges(6,2) + t124 * t4 + t99) + g(2) * (t1 * pkin(9) + (t1 * t57 - t52 * t60) * rSges(6,1) + (-t1 * t52 - t57 * t60) * rSges(6,2) + t124 * t2 + t100) + (t101 + (t57 * rSges(6,1) - t52 * rSges(6,2) + pkin(9)) * t5 + t124 * t6) * g(3)) - m(7) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(7,1) * t24 + rSges(7,2) * t23 + pkin(14) + t112) + (-g(2) * rSges(7,3) + g(1) * t93) * t60 + (g(1) * rSges(7,3) + g(2) * t93) * t55) - m(8) * (g(1) * (rSges(8,1) * t17 - rSges(8,2) * t18 + rSges(8,3) * t55 + t104) + g(2) * (rSges(8,1) * t15 - rSges(8,2) * t16 - rSges(8,3) * t60 + t105) + g(3) * (rSges(8,1) * t19 + rSges(8,2) * t20 + t106)) - m(9) * (g(1) * ((t10 * t17 - t18 * t9) * rSges(9,1) + (-t10 * t18 - t17 * t9) * rSges(9,2) + t55 * rSges(9,3) + t104) + g(2) * ((t10 * t15 - t16 * t9) * rSges(9,1) + (-t10 * t16 - t15 * t9) * rSges(9,2) - t60 * rSges(9,3) + t105) + g(3) * ((t10 * t19 + t20 * t9) * rSges(9,1) + (t10 * t20 - t19 * t9) * rSges(9,2) + t106) + (g(1) * (t17 * t51 + t18 * t50) + g(2) * (t15 * t51 + t16 * t50) + g(3) * (t19 * t51 - t20 * t50)) * pkin(3));
U = t11;
