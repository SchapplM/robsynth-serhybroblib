% Calculate Gravitation load on the joints for
% palh1m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 20:03
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 20:02:26
% EndTime: 2020-04-15 20:02:28
% DurationCPUTime: 1.73s
% Computational Cost: add. (1004->205), mult. (901->287), div. (20->8), fcn. (707->38), ass. (0->124)
t95 = sin(qJ(5));
t155 = rSges(6,2) * t95;
t184 = pkin(11) + rSges(6,3);
t139 = qJ(3) + qJ(4);
t86 = qJ(2) + t139;
t68 = sin(t86);
t70 = cos(t86);
t185 = t68 * t155 + t70 * t184;
t88 = -qJ(7) + pkin(19);
t73 = -qJ(10) + t88;
t64 = -qJ(2) + t73;
t48 = sin(t64);
t49 = cos(t64);
t183 = rSges(11,1) * t49 + rSges(11,2) * t48;
t90 = qJ(8) + qJ(9);
t85 = qJ(2) + t90;
t67 = sin(t85);
t69 = cos(t85);
t146 = t67 * rSges(10,1) + t69 * rSges(10,2);
t91 = qJ(2) + qJ(8);
t78 = sin(t91);
t182 = -pkin(2) * t78 + t146;
t172 = t70 * pkin(9) + t184 * t68;
t101 = cos(qJ(1));
t97 = sin(qJ(1));
t174 = g(1) * t101 + g(2) * t97;
t99 = cos(qJ(5));
t160 = rSges(6,1) * t99;
t178 = (-pkin(9) - t160) * t68;
t102 = t174 * t178;
t136 = t185 * t101;
t137 = t185 * t97;
t103 = g(1) * t136 + g(2) * t137;
t109 = (-t155 + t160) * t70 + t172;
t123 = -rSges(5,1) * t68 - rSges(5,2) * t70;
t110 = t123 * t101;
t111 = t123 * t97;
t177 = t70 * rSges(5,1) - rSges(5,2) * t68;
t180 = pkin(8) * (-m(5) * (g(1) * t110 + g(2) * t111 + g(3) * t177) - m(6) * (g(3) * t109 + t102 + t103));
t179 = g(2) * t101;
t132 = -rSges(11,1) * t48 + t49 * rSges(11,2);
t76 = -qJ(2) + t88;
t129 = t132 + pkin(4) * sin(t76);
t93 = qJ(2) + qJ(3);
t80 = sin(t93);
t84 = cos(t93);
t176 = t84 * rSges(4,1) - rSges(4,2) * t80;
t147 = t183 * t101;
t148 = t183 * t97;
t173 = pkin(10) * m(11) * (g(1) * t147 + g(2) * t148 + g(3) * t132);
t96 = sin(qJ(2));
t169 = pkin(1) * t96;
t167 = pkin(4) * cos(t76);
t166 = pkin(5) * t80;
t66 = pkin(15) - t169;
t165 = t66 * t179;
t164 = pkin(15) * t179;
t100 = cos(qJ(2));
t161 = pkin(1) * t100;
t92 = qJ(2) + qJ(7);
t83 = cos(t92);
t159 = rSges(8,1) * t83;
t158 = rSges(10,1) * t69;
t79 = sin(t92);
t154 = rSges(8,2) * t79;
t153 = rSges(10,2) * t67;
t150 = t97 * t95;
t149 = t97 * t99;
t141 = t101 * t95;
t140 = t101 * t99;
t82 = cos(t91);
t118 = -rSges(9,1) * t78 - rSges(9,2) * t82;
t128 = -pkin(2) * t82 - t153;
t35 = t97 * t158;
t38 = t101 * t158;
t131 = -m(9) * (g(3) * t118 + t174 * (-rSges(9,1) * t82 + rSges(9,2) * t78)) - m(10) * (g(1) * (t128 * t101 + t38) + g(2) * (t128 * t97 + t35) + g(3) * t182);
t65 = pkin(5) * t84;
t130 = t177 + t65;
t127 = -t159 - t161;
t125 = -rSges(4,1) * t80 - rSges(4,2) * t84;
t94 = sin(qJ(6));
t98 = cos(qJ(6));
t122 = rSges(7,1) * t98 - rSges(7,2) * t94;
t120 = rSges(8,1) * t79 + rSges(8,2) * t83;
t116 = -rSges(3,1) * t96 - rSges(3,2) * t100;
t114 = -pkin(14) + t122;
t113 = pkin(15) + t182;
t112 = -t66 - t129;
t104 = t109 + t65;
t89 = qJ(3) + pkin(17);
t81 = cos(t90);
t77 = sin(t90);
t75 = pkin(20) + t92;
t74 = pkin(18) + t139;
t72 = cos(t89);
t71 = sin(t89);
t62 = cos(t75);
t61 = cos(t74);
t60 = sin(t75);
t59 = sin(t74);
t57 = cos(t73);
t56 = sin(t73);
t45 = t101 * t154;
t44 = t97 * t154;
t34 = -t161 - t166;
t33 = t81 * pkin(12) - cos(qJ(8)) * pkin(2);
t32 = -t77 * pkin(12) + sin(qJ(8)) * pkin(2);
t27 = t65 + t66;
t26 = pkin(3) * t62 + t161;
t25 = pkin(3) * t60 + t169;
t24 = -t161 - t167;
t23 = pkin(10) * t61 + cos(qJ(3)) * pkin(5);
t22 = pkin(10) * t59 + sin(qJ(3)) * pkin(5);
t18 = t101 * t34;
t17 = t97 * t34;
t16 = t57 * pkin(8) - pkin(4) * cos(t88);
t15 = t56 * pkin(8) - pkin(4) * sin(t88);
t14 = t101 * t27;
t13 = t70 * t140 + t150;
t12 = -t70 * t141 + t149;
t11 = -t70 * t149 + t141;
t10 = t70 * t150 + t140;
t8 = 0.1e1 / (-t56 * t59 + t57 * t61) / pkin(10) / pkin(8);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t97 - rSges(2,2) * t101) + g(2) * (rSges(2,1) * t101 - rSges(2,2) * t97)) - m(3) * (t164 + (g(1) * rSges(3,3) + g(2) * t116) * t101 + (g(1) * (-pkin(15) - t116) + g(2) * rSges(3,3)) * t97) - m(4) * (t165 + (g(1) * rSges(4,3) + g(2) * t176) * t101 + (g(1) * (-t176 - t66) + g(2) * rSges(4,3)) * t97) - m(5) * (g(2) * t14 + (g(1) * rSges(5,3) + g(2) * t177) * t101 + (g(1) * (-t177 - t27) + g(2) * rSges(5,3)) * t97) - m(6) * ((t13 * rSges(6,1) + t12 * rSges(6,2) + t101 * t172 + t14) * g(2) + (rSges(6,1) * t11 + rSges(6,2) * t10 + (-t172 - t27) * t97) * g(1)) - m(7) * ((g(2) * rSges(7,3) - g(1) * t114) * t97 + (g(1) * rSges(7,3) + g(2) * t114) * t101) - m(8) * (t165 + (g(1) * rSges(8,3) - g(2) * t120) * t101 + (g(1) * (t120 - t66) + g(2) * rSges(8,3)) * t97) - m(9) * (t164 + (g(1) * rSges(9,3) + g(2) * t118) * t101 + (g(1) * (-pkin(15) - t118) + g(2) * rSges(9,3)) * t97) - m(10) * ((g(2) * rSges(10,3) - g(1) * t113) * t97 + (g(1) * rSges(10,3) + g(2) * t113) * t101) - m(11) * ((g(2) * rSges(11,3) + g(1) * t112) * t97 + (g(1) * rSges(11,3) - g(2) * t112) * t101); -m(3) * (g(3) * t116 + t174 * (-rSges(3,1) * t100 + rSges(3,2) * t96)) - m(4) * (g(3) * (t176 - t169) + t174 * (t125 - t161)) - m(5) * (g(1) * (t18 + t110) + g(2) * (t17 + t111) + g(3) * (t130 - t169)) - m(6) * (g(1) * (t18 + t136) + g(2) * (t17 + t137) + g(3) * (t104 - t169) + t102) - m(8) * (g(1) * (t101 * t127 + t45) + g(2) * (t127 * t97 + t44) + g(3) * (-t120 - t169)) - m(11) * (g(1) * (t101 * t24 + t147) + g(2) * (t24 * t97 + t148) + g(3) * (t129 - t169)) + t131 + (-(t25 * t62 - t26 * t60) * pkin(3) * m(7) * (g(3) * t122 + t174 * (-rSges(7,1) * t94 - rSges(7,2) * t98)) + (-m(8) * (g(1) * (-t101 * t159 + t45) + g(2) * (-t159 * t97 + t44) - g(3) * t120) - m(11) * (g(1) * (-t101 * t167 + t147) + g(2) * (-t167 * t97 + t148) + g(3) * t129) + (-(t15 * t57 - t16 * t56) * t180 + (-t15 * t59 + t16 * t61) * t173) * t8) * (t25 * t94 + t26 * t98) * pkin(7)) / (-t60 * t94 - t62 * t98) / pkin(7) / pkin(3); -m(4) * (g(3) * t176 + t125 * t174) - m(5) * (g(3) * t130 + t174 * (t123 - t166)) - m(6) * (g(3) * t104 + t103 + t174 * (-t166 + t178)) + ((t22 * t56 - t23 * t57) * t180 - (-t22 * t61 + t23 * t59) * t173) * t8 + ((-t71 * t77 - t72 * t81) * pkin(12) * t131 - (-t32 * t71 + t33 * t72) * m(10) * (g(1) * (-t101 * t153 + t38) + g(2) * (-t153 * t97 + t35) + g(3) * t146)) * pkin(6) / (t32 * t81 + t33 * t77) / pkin(12); -m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t13) + g(2) * (-rSges(6,1) * t10 + rSges(6,2) * t11) + g(3) * (-rSges(6,1) * t95 - rSges(6,2) * t99) * t68);];
taug = t1(:);
