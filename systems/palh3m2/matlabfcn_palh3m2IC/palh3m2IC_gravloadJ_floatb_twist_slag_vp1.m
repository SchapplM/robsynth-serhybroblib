% Calculate Gravitation load on the joints for
% palh3m2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:59:43
% EndTime: 2020-05-07 04:59:45
% DurationCPUTime: 1.45s
% Computational Cost: add. (768->198), mult. (607->234), div. (14->8), fcn. (492->39), ass. (0->101)
t135 = pkin(10) + rSges(6,3);
t58 = qJ(2) + qJ(3);
t55 = qJ(4) + t58;
t47 = sin(t55);
t154 = t135 * t47;
t60 = sin(qJ(5));
t128 = rSges(6,2) * t60;
t64 = cos(qJ(5));
t132 = rSges(6,1) * t64;
t153 = (-t128 + t132) * t47;
t112 = (qJ(2) - qJ(6));
t108 = (qJ(4) + pkin(14));
t100 = (qJ(3) + t108);
t98 = (pkin(15) + t100);
t78 = (t98 - t112);
t72 = -2 * qJ(7) - pkin(16) + t78;
t79 = (t98 + t112);
t73 = pkin(16) + t79;
t152 = -cos((qJ(8) - t72)) + cos((qJ(8) - t73));
t111 = pkin(15) - qJ(7);
t54 = -qJ(2) + t111;
t105 = -qJ(8) + t111;
t49 = -qJ(2) + t105;
t41 = sin(t49);
t42 = cos(t49);
t88 = rSges(9,1) * t42 + rSges(9,2) * t41;
t82 = pkin(3) * cos(t54) - t88;
t48 = cos(t55);
t151 = -rSges(5,1) * t48 + t47 * rSges(5,2);
t51 = sin(t58);
t53 = cos(t58);
t150 = -rSges(4,1) * t53 + t51 * rSges(4,2);
t57 = qJ(2) + qJ(7);
t50 = sin(t57);
t52 = cos(t57);
t149 = t52 * rSges(8,1) - rSges(8,2) * t50;
t131 = rSges(9,1) * t41;
t97 = -t131 + pkin(3) * sin(t54);
t62 = sin(qJ(1));
t66 = cos(qJ(1));
t148 = -g(1) * t66 - g(2) * t62;
t147 = -t48 * pkin(8) - t154;
t61 = sin(qJ(2));
t144 = pkin(1) * t61;
t142 = pkin(4) * t53;
t65 = cos(qJ(2));
t56 = t65 * pkin(1);
t46 = t56 + pkin(12);
t140 = g(2) * t66 * t46;
t124 = t47 * t66;
t129 = rSges(5,2) * t48;
t116 = rSges(5,1) * t124 + t66 * t129;
t125 = t47 * t62;
t117 = rSges(5,1) * t125 + t62 * t129;
t71 = (g(3) * (-pkin(8) - t132) + t148 * t135) * t48;
t80 = pkin(8) * t124 + t153 * t66;
t81 = pkin(8) * t125 + t153 * t62;
t87 = t48 * t128 - t154;
t138 = (-m(5) * (g(1) * t116 + g(2) * t117 + g(3) * t151) - m(6) * (g(1) * t80 + g(2) * t81 + g(3) * t87 + t71)) / pkin(9);
t130 = rSges(4,2) * t53;
t126 = rSges(9,2) * t42;
t123 = t51 * t62;
t122 = t51 * t66;
t121 = t60 * t62;
t120 = t60 * t66;
t119 = t62 * t64;
t118 = t64 * t66;
t115 = rSges(4,1) * t123 + t62 * t130;
t114 = rSges(4,1) * t122 + t66 * t130;
t110 = qJ(7) + pkin(16);
t18 = t62 * t126;
t19 = t66 * t126;
t109 = m(9) * (g(1) * (-t66 * t131 + t19) + g(2) * (-t62 * t131 + t18) - g(3) * t88) / pkin(7);
t102 = -t144 + t97;
t99 = t110 + t112;
t96 = rSges(3,1) * t65 - rSges(3,2) * t61;
t59 = sin(qJ(6));
t63 = cos(qJ(6));
t92 = rSges(7,1) * t63 - rSges(7,2) * t59;
t89 = -rSges(8,1) * t50 - rSges(8,2) * t52;
t86 = t151 - t142;
t85 = -pkin(6) + t92;
t84 = pkin(12) + t96;
t83 = -t46 - t82;
t76 = -qJ(7) + t78;
t75 = -qJ(7) + t79;
t74 = t87 - t142;
t70 = 0.1e1 / pkin(2);
t40 = sin(t99);
t35 = pkin(4) * t122;
t34 = pkin(4) * t123;
t20 = pkin(4) * t51 - t144;
t17 = t46 - t142;
t10 = t66 * t20;
t9 = t62 * t20;
t8 = t66 * t17;
t7 = t48 * t118 - t121;
t6 = t48 * t120 + t119;
t5 = t48 * t119 + t120;
t4 = t48 * t121 - t118;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t62 - rSges(2,2) * t66) + g(2) * (rSges(2,1) * t66 - rSges(2,2) * t62)) - m(3) * ((g(1) * rSges(3,3) + g(2) * t84) * t66 + (g(2) * rSges(3,3) - g(1) * t84) * t62) - m(4) * (t140 + (g(1) * rSges(4,3) + g(2) * t150) * t66 + (g(1) * (-t46 - t150) + g(2) * rSges(4,3)) * t62) - m(5) * (g(2) * t8 + (g(1) * rSges(5,3) + g(2) * t151) * t66 + (g(1) * (-t17 - t151) + g(2) * rSges(5,3)) * t62) - m(6) * ((-rSges(6,1) * t7 + rSges(6,2) * t6 + t147 * t66 + t8) * g(2) + (rSges(6,1) * t5 - rSges(6,2) * t4 + (-t147 - t17) * t62) * g(1)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t85) * t66 + (g(2) * rSges(7,3) - g(1) * t85) * t62) - m(8) * (t140 + (g(1) * rSges(8,3) + g(2) * t149) * t66 + (g(1) * (-t46 - t149) + g(2) * rSges(8,3)) * t62) - m(9) * ((g(1) * rSges(9,3) - g(2) * t83) * t66 + (g(2) * rSges(9,3) + g(1) * t83) * t62); -m(3) * (g(3) * t96 - t148 * (-rSges(3,1) * t61 - rSges(3,2) * t65)) - m(4) * (g(1) * (-t144 * t66 + t114) + g(2) * (-t144 * t62 + t115) + g(3) * (t150 + t56)) - m(5) * (g(1) * (t10 + t116) + g(2) * (t9 + t117) + g(3) * (t56 + t86)) - m(6) * (g(1) * (t10 + t80) + g(2) * (t81 + t9) + g(3) * (t56 + t74) + t71) - m(8) * (g(3) * (t149 + t56) - t148 * (t89 - t144)) - m(9) * (g(1) * (t102 * t66 + t19) + g(2) * (t102 * t62 + t18) + g(3) * (t56 + t82)) + (-pkin(1) * sin(t110) / pkin(5) * m(7) * (g(3) * t92 - t148 * (-rSges(7,1) * t59 - rSges(7,2) * t63)) + (-pkin(2) * t40 - pkin(1) * sin(t112)) * t70 * (-m(8) * (g(3) * t149 - t148 * t89) - m(9) * (g(1) * (t97 * t66 + t19) + g(2) * (t97 * t62 + t18) + g(3) * t82))) / t40 + (pkin(3) * (pkin(1) * (cos((-qJ(8) + t112)) - cos((qJ(8) + t112))) + (cos((-qJ(8) + t99)) - cos((qJ(8) + t99))) * pkin(2)) * t138 - (((cos(t76) - cos(t75)) * pkin(1) + (cos(t72) - cos(t73)) * pkin(2)) * pkin(3) + ((-cos((-qJ(8) + t76)) + cos((-qJ(8) + t75))) * pkin(1) + t152 * pkin(2)) * pkin(7)) * t109) / t152 * t70; -m(4) * (g(1) * t114 + g(2) * t115 + g(3) * t150) - m(5) * (g(1) * (t35 + t116) + g(2) * (t34 + t117) + g(3) * t86) - m(6) * (g(1) * (t35 + t80) + g(2) * (t34 + t81) + g(3) * t74 + t71) - pkin(9) * t138 + (-(cos((2 * qJ(3) + t108)) - cos((2 * qJ(7) - 2 * pkin(15) + 2 * qJ(8) + t108))) / (cos((2 * t100)) - cos((2 * t105))) * t138 - sin(t108) / sin((-qJ(7) - qJ(8) + t98)) * t109) * pkin(4); -m(6) * (g(1) * (rSges(6,1) * t6 + rSges(6,2) * t7) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5) + g(3) * (rSges(6,1) * t60 + rSges(6,2) * t64) * t47);];
taug = t1(:);
