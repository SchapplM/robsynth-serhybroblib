% Calculate Gravitation load on the joints for
% picker2Dm2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% rSges [11x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm2IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:20:52
% EndTime: 2020-05-11 09:20:53
% DurationCPUTime: 0.93s
% Computational Cost: add. (998->178), mult. (466->188), div. (37->10), fcn. (360->59), ass. (0->113)
t76 = qJ(1) + qJ(2);
t64 = sin(t76);
t66 = cos(t76);
t77 = sin(qJ(7));
t81 = cos(qJ(7));
t143 = (t64 * t77 + t66 * t81) * pkin(3);
t69 = qJ(6) + t76;
t50 = sin(t69);
t53 = cos(t69);
t104 = -rSges(7,1) * t50 - rSges(7,2) * t53;
t114 = t53 * rSges(7,1) - rSges(7,2) * t50;
t142 = m(7) * (g(1) * t104 + g(2) * t114);
t121 = qJ(9) + t76;
t59 = qJ(3) + t121;
t38 = sin(t59);
t41 = cos(t59);
t116 = -rSges(10,1) * t41 + t38 * rSges(10,2);
t131 = t38 * rSges(10,1) + t41 * rSges(10,2);
t141 = m(10) * (g(1) * t131 + g(2) * t116);
t82 = cos(qJ(1));
t140 = pkin(1) * t82;
t139 = pkin(2) * t66;
t138 = pkin(3) * t66;
t125 = qJ(2) + qJ(4);
t70 = qJ(1) + t125;
t54 = cos(t70);
t137 = pkin(4) * t54;
t56 = qJ(10) + t70;
t33 = sin(t56);
t34 = cos(t56);
t112 = t34 * rSges(11,1) - rSges(11,2) * t33;
t98 = -rSges(11,1) * t33 - rSges(11,2) * t34;
t136 = m(11) * (g(1) * t98 + g(2) * t112);
t83 = 0.1e1 / pkin(5);
t124 = qJ(3) - qJ(6);
t73 = pkin(8) + qJ(5);
t106 = t73 + t124;
t75 = qJ(1) + qJ(8);
t89 = t106 - t75;
t107 = t73 - t124;
t90 = t107 - t75;
t135 = 0.1e1 / (pkin(6) * (cos(t90) - cos(t89)) + (-cos(qJ(9) - t90) + cos(qJ(9) + t89)) * pkin(2)) * t83;
t120 = -qJ(9) - t124;
t8 = 0.1e1 / ((-cos(qJ(4) - t124) + cos(qJ(4) + t124)) * pkin(6) + (cos(qJ(4) + t120) - cos(qJ(4) - t120)) * pkin(2));
t85 = 0.1e1 / pkin(3);
t134 = t8 * t85;
t58 = -qJ(9) + t70;
t23 = t58 - t124;
t57 = qJ(4) + t121;
t24 = t57 + t124;
t133 = sin(t24) - sin(t23);
t132 = cos(t24) - cos(t23);
t51 = sin(t70);
t130 = t51 * rSges(5,1) + t54 * rSges(5,2);
t129 = -cos(0.2e1 * t75) + cos(0.2e1 * t73);
t128 = sin(t58) - sin(t57);
t127 = cos(t58) - cos(t57);
t126 = t64 * rSges(3,1) + t66 * rSges(3,2);
t78 = sin(qJ(4));
t74 = 0.1e1 / t78;
t84 = 0.1e1 / pkin(4);
t123 = pkin(1) * t74 * t84;
t48 = pkin(3) * t64;
t122 = t48 + t130;
t119 = -qJ(1) + t73;
t118 = -rSges(3,1) * t66 + t64 * rSges(3,2);
t117 = -rSges(5,1) * t54 + t51 * rSges(5,2);
t71 = qJ(3) + t76;
t52 = sin(t71);
t55 = cos(t71);
t115 = t55 * rSges(4,1) - rSges(4,2) * t52;
t63 = sin(t75);
t65 = cos(t75);
t113 = t65 * rSges(9,1) - rSges(9,2) * t63;
t111 = -pkin(6) * t52 + t131;
t110 = pkin(6) * t55 + t116;
t109 = -qJ(4) + t119;
t108 = qJ(4) + t119;
t105 = -rSges(4,1) * t52 - rSges(4,2) * t55;
t103 = -rSges(9,1) * t63 - rSges(9,2) * t65;
t102 = t51 * t66 - t54 * t64;
t101 = t51 * t77 + t54 * t81;
t49 = pkin(2) * t64;
t99 = t111 + t49;
t97 = -qJ(8) + t109;
t96 = -qJ(8) + t108;
t95 = t115 - t139;
t94 = t117 - t138;
t93 = t112 - t137;
t92 = t105 + t49;
t32 = pkin(4) * t51;
t91 = t32 + t98;
t88 = t110 - t139;
t87 = t48 + t91;
t86 = t93 - t138;
t80 = sin(qJ(1));
t79 = sin(qJ(2));
t72 = t80 * pkin(1);
t68 = -qJ(9) + t73;
t67 = qJ(9) + t73;
t62 = sin(t125);
t61 = cos(t73);
t60 = sin(t73);
t43 = qJ(9) + t106;
t42 = -qJ(9) + t107;
t14 = -pkin(5) * t63 + t72;
t13 = -pkin(5) * t65 + t140;
t12 = t137 + t138 + t140;
t11 = t72 + t48 + t32;
t3 = -m(4) * (g(1) * t105 + g(2) * t115) - m(10) * (g(1) * t111 + g(2) * t110);
t2 = -m(5) * (g(1) * t130 + g(2) * t117) - m(11) * (g(1) * t91 + g(2) * t93);
t1 = -m(3) * (g(1) * t126 + g(2) * t118) - m(4) * (g(1) * t92 + g(2) * t95) - m(5) * (g(1) * t122 + g(2) * t94) - t142 - m(10) * (g(1) * t99 + g(2) * t88) - m(11) * (g(1) * t87 + g(2) * t86);
t4 = [-m(2) * (g(1) * (rSges(2,1) * t80 + rSges(2,2) * t82) + g(2) * (-rSges(2,1) * t82 + rSges(2,2) * t80)) - m(3) * (g(1) * (t72 + t126) + g(2) * (t118 - t140)) - m(4) * (g(1) * (t72 + t92) + g(2) * (t95 - t140)) - m(5) * (g(1) * (t72 + t122) + g(2) * (t94 - t140)) - m(7) * (g(1) * (t104 + t72) + g(2) * (t114 - t140)) - m(10) * (g(1) * (t72 + t99) + g(2) * (t88 - t140)) - m(11) * (g(1) * (t72 + t87) + g(2) * (t86 - t140)) - sin(qJ(8)) / sin(-qJ(8) + t119) * m(6) * (g(1) * (-rSges(6,1) * t60 - rSges(6,2) * t61) + g(2) * (rSges(6,1) * t61 - rSges(6,2) * t60)) + t79 * t123 * t136 + (-g(1) * (t103 + t72) - g(2) * (t113 - t140) - (-t129 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) * t83 / t129 * (g(1) * t103 + g(2) * t113)) * m(9) + (((-t133 * t11 - t132 * t12) * t134 + ((sin(t43) - sin(t42)) * t14 + (cos(t43) - cos(t42)) * t13) * t135) * t3 - ((t128 * t11 + t127 * t12) * t134 + (-(sin(t68) - sin(t67)) * t14 - (cos(t68) - cos(t67)) * t13) * t135) * t142) * pkin(2) + ((-pkin(1) * t62 - pkin(3) * t78) * t74 * t1 + (pkin(3) * t79 + pkin(4) * t62) * t2 * t123 + (pkin(3) * (cos(t109) - cos(t108)) + (cos(-qJ(2) + t97) - cos(qJ(2) + t96)) * pkin(5)) * pkin(1) * t83 / (cos(t97) - cos(t96)) * t141) * t85; -m(8) * (g(1) * (rSges(8,1) * t81 - rSges(8,2) * t77) + g(2) * (rSges(8,1) * t77 + rSges(8,2) * t81)) + ((-t132 * t77 + t133 * t81) * t3 - (t127 * t77 - t128 * t81) * t142) * t8 * pkin(2) + (((t101 * pkin(4) + t143) * t2 - (t102 * pkin(4) - t143) * t136) * t84 - (t1 - t141) * t101) / t102;];
taug = t4(:);
