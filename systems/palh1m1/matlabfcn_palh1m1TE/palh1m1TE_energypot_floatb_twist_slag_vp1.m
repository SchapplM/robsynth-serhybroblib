% Calculate potential energy for
% palh1m1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [23x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DA,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1TE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:10:21
% EndTime: 2020-04-12 20:10:28
% DurationCPUTime: 6.04s
% Computational Cost: add. (80084->239), mult. (120425->345), div. (5724->9), fcn. (76231->28), ass. (0->141)
t154 = -pkin(12) - rSges(6,3);
t82 = pkin(6) ^ 2;
t85 = pkin(1) ^ 2;
t146 = t82 + t85;
t65 = sin(pkin(20));
t67 = cos(pkin(20));
t69 = sin(qJ(3));
t74 = cos(qJ(3));
t157 = pkin(6) * (-t65 * t74 - t67 * t69);
t139 = pkin(1) * t157;
t51 = -0.2e1 * t139;
t84 = 0.1e1 / pkin(2);
t152 = 0.1e1 / (t51 + t146) * t84;
t70 = sin(qJ(2));
t165 = -t70 / 0.2e1;
t156 = pkin(6) * (t65 * t69 - t67 * t74);
t148 = t51 + t82;
t159 = pkin(13) - pkin(2);
t164 = -pkin(2) - pkin(13);
t37 = sqrt(-((pkin(1) - t159) * (pkin(1) + t159) + t148) * ((pkin(1) - t164) * (pkin(1) + t164) + t148));
t131 = -pkin(13) ^ 2 + t146;
t83 = pkin(2) ^ 2;
t40 = t51 + t83 + t131;
t50 = -pkin(1) + t157;
t35 = -t156 * t37 - t40 * t50;
t36 = t156 * t40 - t37 * t50;
t75 = cos(qJ(2));
t170 = (t75 * t35 / 0.2e1 + t36 * t165) * t152;
t169 = 0.1e1 / pkin(3);
t168 = -t37 / 0.2e1;
t167 = t37 / 0.2e1;
t166 = t131 / 0.2e1 - t83 / 0.2e1 - t139;
t163 = pkin(3) - pkin(8);
t162 = -pkin(8) - pkin(3);
t161 = -pkin(9) - pkin(11);
t160 = -pkin(9) + pkin(11);
t158 = pkin(1) * t70;
t72 = sin(pkin(19));
t77 = cos(pkin(19));
t57 = t70 * t77 - t72 * t75;
t155 = pkin(7) * t57;
t153 = sin(pkin(18));
t147 = -0.2e1 * pkin(1) * t155 + t85;
t132 = pkin(7) ^ 2 + t147;
t43 = 0.1e1 / t132;
t151 = t43 / pkin(8);
t38 = sqrt(-((pkin(7) - t163) * (pkin(7) + t163) + t147) * ((pkin(7) - t162) * (pkin(7) + t162) + t147));
t60 = t70 * t72 + t75 * t77;
t150 = t60 * t38;
t149 = 0.1e1 / pkin(13) * t84;
t145 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t144 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t143 = cos(pkin(21));
t142 = cos(pkin(23));
t141 = sin(pkin(21));
t140 = sin(pkin(23));
t138 = pkin(14) + r_base(3);
t71 = sin(qJ(1));
t137 = t71 * pkin(16) + r_base(2);
t76 = cos(qJ(1));
t136 = t76 * pkin(16) + r_base(1);
t135 = rSges(10,1) * t149;
t134 = rSges(10,2) * t149;
t133 = pkin(23) + pkin(22);
t130 = cos(pkin(18)) / 0.2e1;
t129 = pkin(1) - t155;
t128 = t75 * pkin(1) + t138;
t127 = cos(t133);
t126 = sin(t133);
t58 = t69 * t75 + t70 * t74;
t125 = t58 * pkin(5) + t128;
t124 = t132 - t145;
t123 = -rSges(3,1) * t70 - rSges(3,2) * t75;
t112 = t169 * (-pkin(7) * t150 + t124 * t129);
t110 = -t112 / 0.2e1;
t113 = t169 * (pkin(7) * t124 * t60 + t129 * t38);
t111 = t113 / 0.2e1;
t27 = (t110 * t142 + t111 * t140) * t43;
t28 = (t142 * t111 + t140 * t112 / 0.2e1) * t43;
t19 = t27 * t75 - t28 * t70;
t20 = -t27 * t70 - t28 * t75;
t59 = -t69 * t70 + t74 * t75;
t41 = t132 + t145;
t52 = pkin(1) * t57 - pkin(7);
t122 = pkin(1) * t41 * t60 - t38 * t52;
t121 = -pkin(1) * t150 - t41 * t52;
t120 = -t158 * t71 + t137;
t119 = -t158 * t76 + t136;
t29 = (t121 * t130 - t153 * t122 / 0.2e1) * t151;
t30 = (t122 * t130 + t121 * t153 / 0.2e1) * t151;
t118 = rSges(7,1) * t29 - rSges(7,2) * t30 - pkin(15);
t46 = t59 * t71;
t116 = t46 * pkin(5) + t120;
t48 = t59 * t76;
t115 = t48 * pkin(5) + t119;
t26 = (-t75 * t36 / 0.2e1 + t35 * t165) * t152;
t109 = t69 * t110 - t74 * t113 / 0.2e1;
t108 = t110 * t74 + t111 * t69;
t107 = t43 * (t108 * t126 + t109 * t127);
t106 = t43 * (t108 * t127 - t109 * t126);
t105 = pkin(5) * t107;
t104 = pkin(4) - t105;
t103 = -pkin(4) * t107 + pkin(5);
t102 = -0.2e1 * pkin(4) * t105 + pkin(5) ^ 2;
t101 = pkin(4) ^ 2 + t102;
t100 = 0.1e1 / t101;
t99 = 0.1e1 / pkin(9) * t100;
t98 = 0.1e1 / pkin(11) * t100;
t97 = t101 + t144;
t96 = t101 - t144;
t95 = sqrt(-((pkin(4) - t160) * (pkin(4) + t160) + t102) * ((pkin(4) - t161) * (pkin(4) + t161) + t102));
t94 = t95 * t106;
t93 = (-pkin(5) * t94 + t104 * t96) * t99;
t92 = (-pkin(4) * t94 + t103 * t97) * t98;
t91 = -(pkin(5) * t106 * t96 + t104 * t95) * t99 / 0.2e1;
t90 = (pkin(4) * t106 * t97 + t103 * t95) * t98 / 0.2e1;
t73 = cos(qJ(4));
t68 = sin(qJ(4));
t66 = cos(pkin(22));
t64 = sin(pkin(22));
t49 = t58 * t76;
t47 = t58 * t71;
t24 = t76 * t170;
t23 = t76 * t26;
t22 = t71 * t170;
t21 = t71 * t26;
t18 = t19 * t76;
t17 = t20 * t76;
t16 = t19 * t71;
t15 = t20 * t71;
t10 = -t143 * t92 / 0.2e1 + t141 * t90;
t9 = t141 * t92 / 0.2e1 + t143 * t90;
t8 = -t66 * t93 / 0.2e1 + t64 * t91;
t7 = t64 * t93 / 0.2e1 + t66 * t91;
t6 = t10 * t58 + t59 * t9;
t5 = t10 * t59 - t58 * t9;
t4 = t10 * t48 - t49 * t9;
t3 = -t10 * t49 - t48 * t9;
t2 = t10 * t46 - t47 * t9;
t1 = -t10 * t47 - t46 * t9;
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t76 - rSges(2,2) * t71 + r_base(1)) + g(2) * (rSges(2,1) * t71 + rSges(2,2) * t76 + r_base(2)) + g(3) * (rSges(2,3) + t138)) - m(3) * (g(1) * (rSges(3,3) * t71 + t123 * t76 + t136) + g(2) * (-rSges(3,3) * t76 + t123 * t71 + t137) + g(3) * (rSges(3,1) * t75 - rSges(3,2) * t70 + t138)) - m(4) * (g(1) * (rSges(4,1) * t48 - rSges(4,2) * t49 + rSges(4,3) * t71 + t119) + g(2) * (rSges(4,1) * t46 - rSges(4,2) * t47 - rSges(4,3) * t76 + t120) + g(3) * (rSges(4,1) * t58 + rSges(4,2) * t59 + t128)) - m(5) * (g(1) * (rSges(5,1) * t4 + rSges(5,2) * t3 + rSges(5,3) * t71 + t115) + g(2) * (rSges(5,1) * t2 + rSges(5,2) * t1 - rSges(5,3) * t76 + t116) + g(3) * (rSges(5,1) * t6 + rSges(5,2) * t5 + t125)) - m(6) * (g(1) * (t4 * pkin(10) + (t4 * t73 + t68 * t71) * rSges(6,1) + (-t4 * t68 + t71 * t73) * rSges(6,2) + t154 * t3 + t115) + g(2) * (t2 * pkin(10) + (t2 * t73 - t68 * t76) * rSges(6,1) + (-t2 * t68 - t73 * t76) * rSges(6,2) + t154 * t1 + t116) + (t125 + (t73 * rSges(6,1) - t68 * rSges(6,2) + pkin(10)) * t6 + t154 * t5) * g(3)) - m(7) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(7,1) * t30 + rSges(7,2) * t29 - pkin(17) + t138) + (-g(2) * rSges(7,3) + g(1) * t118) * t76 + (g(1) * rSges(7,3) + g(2) * t118) * t71) - m(8) * (g(1) * (rSges(8,1) * t17 - rSges(8,2) * t18 + rSges(8,3) * t71 + t119) + g(2) * (rSges(8,1) * t15 - rSges(8,2) * t16 - rSges(8,3) * t76 + t120) + g(3) * (rSges(8,1) * t19 + rSges(8,2) * t20 + t128)) - m(9) * (g(1) * (rSges(9,1) * t23 - rSges(9,2) * t24 + rSges(9,3) * t71 + t136) + g(2) * (rSges(9,1) * t21 - rSges(9,2) * t22 - rSges(9,3) * t76 + t137) + g(3) * (rSges(9,1) * t170 + rSges(9,2) * t26 + t138)) - m(10) * (g(1) * (t23 * pkin(2) + (t166 * t23 - t168 * t24) * t135 + (-t166 * t24 + t167 * t23) * t134 + t71 * rSges(10,3) + t136) + g(2) * (t21 * pkin(2) + (t166 * t21 - t168 * t22) * t135 + (-t166 * t22 + t167 * t21) * t134 - t76 * rSges(10,3) + t137) + g(3) * (t170 * pkin(2) + (t166 * t170 + t168 * t26) * t135 + (t166 * t26 + t167 * t170) * t134 + t138)) - m(11) * (g(1) * ((t17 * t8 - t18 * t7) * rSges(11,1) + (-t17 * t7 - t18 * t8) * rSges(11,2) + t71 * rSges(11,3) + t119) + g(2) * ((t15 * t8 - t16 * t7) * rSges(11,1) + (-t15 * t7 - t16 * t8) * rSges(11,2) - t76 * rSges(11,3) + t120) + g(3) * ((t19 * t8 + t20 * t7) * rSges(11,1) + (-t19 * t7 + t20 * t8) * rSges(11,2) + t128) + (g(1) * (t17 * t66 + t18 * t64) + g(2) * (t15 * t66 + t16 * t64) + g(3) * (t19 * t66 - t20 * t64)) * pkin(4));
U = t11;
