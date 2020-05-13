% Calculate Gravitation load on the joints for
% palh3m1IC
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
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m1IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3,1),zeros(16,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [10x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:31:36
% EndTime: 2020-04-20 17:31:37
% DurationCPUTime: 1.36s
% Computational Cost: add. (842->174), mult. (741->242), div. (16->6), fcn. (587->28), ass. (0->102)
t131 = pkin(10) + rSges(6,3);
t109 = qJ(3) + qJ(4);
t67 = qJ(2) + t109;
t56 = sin(t67);
t149 = t131 * t56;
t73 = sin(qJ(5));
t124 = rSges(6,2) * t73;
t77 = cos(qJ(5));
t128 = rSges(6,1) * t77;
t148 = (-t124 + t128) * t56;
t79 = cos(qJ(1));
t120 = t56 * t79;
t57 = cos(t67);
t125 = rSges(5,2) * t57;
t112 = rSges(5,1) * t120 + t79 * t125;
t75 = sin(qJ(1));
t121 = t56 * t75;
t113 = rSges(5,1) * t121 + t75 * t125;
t145 = -rSges(5,1) * t57 + t56 * rSges(5,2);
t142 = -g(1) * t79 - g(2) * t75;
t80 = (g(3) * (-pkin(8) - t128) + t142 * t131) * t57;
t84 = pkin(8) * t120 + t148 * t79;
t85 = pkin(8) * t121 + t148 * t75;
t91 = t57 * t124 - t149;
t147 = pkin(7) * (-m(5) * (g(1) * t112 + g(2) * t113 + g(3) * t145) - m(6) * (g(1) * t84 + g(2) * t85 + g(3) * t91 + t80));
t69 = -qJ(7) + pkin(15);
t66 = -qJ(2) + t69;
t64 = -qJ(8) + t69;
t58 = -qJ(2) + t64;
t44 = sin(t58);
t45 = cos(t58);
t93 = rSges(9,1) * t45 + rSges(9,2) * t44;
t86 = pkin(3) * cos(t66) - t93;
t71 = qJ(2) + qJ(3);
t61 = sin(t71);
t63 = cos(t71);
t144 = -rSges(4,1) * t63 + t61 * rSges(4,2);
t70 = qJ(2) + qJ(7);
t60 = sin(t70);
t62 = cos(t70);
t143 = t62 * rSges(8,1) - rSges(8,2) * t60;
t127 = rSges(9,1) * t44;
t102 = -t127 + pkin(3) * sin(t66);
t122 = rSges(9,2) * t45;
t25 = t75 * t122;
t26 = t79 * t122;
t141 = pkin(9) * m(9) * (g(1) * (-t127 * t79 + t26) + g(2) * (-t127 * t75 + t25) - g(3) * t93);
t140 = -t57 * pkin(8) - t149;
t74 = sin(qJ(2));
t138 = pkin(1) * t74;
t136 = pkin(4) * t63;
t78 = cos(qJ(2));
t68 = t78 * pkin(1);
t55 = t68 + pkin(12);
t134 = g(2) * t79 * t55;
t126 = rSges(4,2) * t63;
t119 = t61 * t75;
t118 = t61 * t79;
t117 = t73 * t79;
t116 = t75 * t73;
t115 = t75 * t77;
t114 = t77 * t79;
t111 = rSges(4,1) * t119 + t75 * t126;
t110 = rSges(4,1) * t118 + t79 * t126;
t104 = -t138 + t102;
t101 = rSges(3,1) * t78 - rSges(3,2) * t74;
t72 = sin(qJ(6));
t76 = cos(qJ(6));
t97 = rSges(7,1) * t76 - rSges(7,2) * t72;
t94 = -rSges(8,1) * t60 - rSges(8,2) * t62;
t90 = t145 - t136;
t89 = -pkin(6) + t97;
t88 = pkin(12) + t101;
t87 = -t55 - t86;
t81 = t91 - t136;
t65 = pkin(16) + t70;
t59 = pkin(14) + t109;
t54 = cos(t65);
t53 = cos(t64);
t51 = sin(t65);
t50 = sin(t64);
t49 = cos(t59);
t48 = sin(t59);
t41 = pkin(4) * t118;
t40 = pkin(4) * t119;
t27 = pkin(4) * t61 - t138;
t24 = t55 - t136;
t23 = pkin(2) * t54 + t68;
t22 = -pkin(2) * t51 - t138;
t20 = pkin(9) * t49 + cos(qJ(3)) * pkin(4);
t19 = -pkin(9) * t48 - sin(qJ(3)) * pkin(4);
t15 = t79 * t27;
t14 = t75 * t27;
t13 = -t53 * pkin(7) + pkin(3) * cos(t69);
t12 = -t50 * pkin(7) + pkin(3) * sin(t69);
t11 = t79 * t24;
t10 = t114 * t57 - t116;
t9 = t117 * t57 + t115;
t8 = t115 * t57 + t117;
t7 = t116 * t57 - t114;
t5 = 0.1e1 / (t48 * t53 + t49 * t50) / pkin(9) / pkin(7);
t1 = [-m(2) * (g(1) * (-t75 * rSges(2,1) - rSges(2,2) * t79) + g(2) * (rSges(2,1) * t79 - t75 * rSges(2,2))) - m(3) * ((g(1) * rSges(3,3) + g(2) * t88) * t79 + (g(2) * rSges(3,3) - g(1) * t88) * t75) - m(4) * (t134 + (g(1) * rSges(4,3) + g(2) * t144) * t79 + (g(1) * (-t55 - t144) + g(2) * rSges(4,3)) * t75) - m(5) * (g(2) * t11 + (g(1) * rSges(5,3) + g(2) * t145) * t79 + (g(1) * (-t24 - t145) + g(2) * rSges(5,3)) * t75) - m(6) * ((-t10 * rSges(6,1) + t9 * rSges(6,2) + t140 * t79 + t11) * g(2) + (rSges(6,1) * t8 - rSges(6,2) * t7 + (-t140 - t24) * t75) * g(1)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t89) * t79 + (g(2) * rSges(7,3) - g(1) * t89) * t75) - m(8) * (t134 + (g(1) * rSges(8,3) + g(2) * t143) * t79 + (g(1) * (-t55 - t143) + g(2) * rSges(8,3)) * t75) - m(9) * ((g(1) * rSges(9,3) - g(2) * t87) * t79 + (g(2) * rSges(9,3) + g(1) * t87) * t75); -m(3) * (g(3) * t101 - t142 * (-rSges(3,1) * t74 - rSges(3,2) * t78)) - m(4) * (g(1) * (-t138 * t79 + t110) + g(2) * (-t138 * t75 + t111) + g(3) * (t144 + t68)) - m(5) * (g(1) * (t15 + t112) + g(2) * (t14 + t113) + g(3) * (t68 + t90)) - m(6) * (g(1) * (t15 + t84) + g(2) * (t14 + t85) + g(3) * (t68 + t81) + t80) - m(8) * (g(3) * (t143 + t68) - t142 * (t94 - t138)) - m(9) * (g(1) * (t104 * t79 + t26) + g(2) * (t104 * t75 + t25) + g(3) * (t68 + t86)) + (-(-t22 * t54 - t23 * t51) * pkin(2) * m(7) * (g(3) * t97 - t142 * (-rSges(7,1) * t72 - rSges(7,2) * t76)) + (m(8) * (g(3) * t143 - t142 * t94) + m(9) * (g(1) * (t102 * t79 + t26) + g(2) * (t102 * t75 + t25) + g(3) * t86) + ((-t12 * t53 + t13 * t50) * t147 - (-t12 * t49 - t13 * t48) * t141) * t5) * pkin(5) * (t22 * t76 + t23 * t72)) / (-t51 * t76 + t54 * t72) / pkin(5) / pkin(2); -m(4) * (g(1) * t110 + g(2) * t111 + g(3) * t144) - m(5) * (g(1) * (t41 + t112) + g(2) * (t40 + t113) + g(3) * t90) - m(6) * (g(1) * (t41 + t84) + g(2) * (t40 + t85) + g(3) * t81 + t80) + ((t19 * t53 - t20 * t50) * t147 - (t19 * t49 + t20 * t48) * t141) * t5; -m(6) * (g(1) * (rSges(6,1) * t9 + rSges(6,2) * t10) + g(2) * (rSges(6,1) * t7 + rSges(6,2) * t8) + g(3) * (rSges(6,1) * t73 + rSges(6,2) * t77) * t56);];
taug = t1(:);
