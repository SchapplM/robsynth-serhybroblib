% Calculate Gravitation load on the joints for
% palh1m1OL
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
% taug [13x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m1OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:36
% EndTime: 2020-04-15 19:28:44
% DurationCPUTime: 1.15s
% Computational Cost: add. (656->174), mult. (604->236), div. (0->0), fcn. (498->22), ass. (0->96)
t64 = sin(qJ(5));
t121 = rSges(6,2) * t64;
t144 = pkin(11) + rSges(6,3);
t62 = qJ(2) + qJ(3);
t58 = qJ(4) + t62;
t47 = sin(t58);
t49 = cos(t58);
t145 = t47 * t121 + t49 * t144;
t61 = qJ(2) + qJ(7);
t50 = pkin(19) - t61;
t43 = -qJ(10) + t50;
t33 = sin(t43);
t34 = cos(t43);
t143 = rSges(11,1) * t34 + rSges(11,2) * t33;
t60 = qJ(2) + qJ(8);
t57 = qJ(9) + t60;
t46 = sin(t57);
t48 = cos(t57);
t108 = t46 * rSges(10,1) + t48 * rSges(10,2);
t51 = sin(t60);
t142 = -pkin(2) * t51 + t108;
t135 = t49 * pkin(9) + t144 * t47;
t70 = cos(qJ(1));
t140 = g(2) * t70;
t98 = -t33 * rSges(11,1) + t34 * rSges(11,2);
t95 = pkin(4) * sin(t50) + t98;
t139 = t49 * rSges(5,1) - t47 * rSges(5,2);
t53 = sin(t62);
t56 = cos(t62);
t138 = t56 * rSges(4,1) - t53 * rSges(4,2);
t66 = sin(qJ(1));
t137 = g(1) * t70 + g(2) * t66;
t68 = cos(qJ(5));
t124 = rSges(6,1) * t68;
t136 = (-pkin(9) - t124) * t47;
t65 = sin(qJ(2));
t134 = pkin(1) * t65;
t69 = cos(qJ(2));
t133 = pkin(1) * t69;
t131 = pkin(4) * cos(t50);
t130 = pkin(5) * t53;
t45 = pkin(15) - t134;
t128 = t45 * t140;
t127 = pkin(15) * t140;
t55 = cos(t61);
t123 = rSges(8,1) * t55;
t122 = rSges(10,1) * t48;
t52 = sin(t61);
t120 = rSges(8,2) * t52;
t119 = rSges(10,2) * t46;
t114 = t64 * t70;
t113 = t66 * t64;
t112 = t66 * t68;
t111 = t68 * t70;
t110 = t143 * t66;
t109 = t143 * t70;
t103 = t145 * t66;
t102 = t145 * t70;
t20 = t66 * t122;
t23 = t70 * t122;
t54 = cos(t60);
t82 = -rSges(9,1) * t51 - rSges(9,2) * t54;
t93 = -pkin(2) * t54 - t119;
t97 = -m(9) * (g(3) * t82 + t137 * (-rSges(9,1) * t54 + rSges(9,2) * t51)) - m(10) * (g(1) * (t70 * t93 + t23) + g(2) * (t66 * t93 + t20) + g(3) * t142);
t44 = pkin(5) * t56;
t96 = t44 + t139;
t94 = -t123 - t133;
t91 = -rSges(3,1) * t65 - rSges(3,2) * t69;
t89 = -rSges(4,1) * t53 - rSges(4,2) * t56;
t87 = -rSges(5,1) * t47 - rSges(5,2) * t49;
t63 = sin(qJ(6));
t67 = cos(qJ(6));
t86 = rSges(7,1) * t67 - rSges(7,2) * t63;
t84 = t52 * rSges(8,1) + t55 * rSges(8,2);
t81 = -pkin(14) + t86;
t80 = pkin(15) + t142;
t79 = -t45 - t95;
t78 = t87 * t66;
t77 = t87 * t70;
t76 = (-t121 + t124) * t49 + t135;
t74 = t44 + t76;
t73 = g(1) * t102 + g(2) * t103;
t72 = t137 * t136;
t30 = t70 * t120;
t29 = t66 * t120;
t19 = -t130 - t133;
t14 = t44 + t45;
t13 = -t131 - t133;
t9 = t70 * t19;
t8 = t66 * t19;
t7 = t70 * t14;
t6 = t111 * t49 + t113;
t5 = -t114 * t49 + t112;
t4 = -t112 * t49 + t114;
t3 = t113 * t49 + t111;
t1 = [-m(2) * (g(1) * (-t66 * rSges(2,1) - rSges(2,2) * t70) + g(2) * (rSges(2,1) * t70 - t66 * rSges(2,2))) - m(3) * (t127 + (g(1) * rSges(3,3) + g(2) * t91) * t70 + (g(1) * (-pkin(15) - t91) + g(2) * rSges(3,3)) * t66) - m(4) * (t128 + (g(1) * rSges(4,3) + g(2) * t138) * t70 + (g(1) * (-t45 - t138) + g(2) * rSges(4,3)) * t66) - m(5) * (g(2) * t7 + (g(1) * rSges(5,3) + g(2) * t139) * t70 + (g(1) * (-t14 - t139) + g(2) * rSges(5,3)) * t66) - m(6) * ((t6 * rSges(6,1) + t5 * rSges(6,2) + t135 * t70 + t7) * g(2) + (t4 * rSges(6,1) + t3 * rSges(6,2) + (-t135 - t14) * t66) * g(1)) - m(7) * ((g(1) * rSges(7,3) + g(2) * t81) * t70 + (g(2) * rSges(7,3) - g(1) * t81) * t66) - m(8) * (t128 + (g(1) * rSges(8,3) - g(2) * t84) * t70 + (g(1) * (-t45 + t84) + g(2) * rSges(8,3)) * t66) - m(9) * (t127 + (g(1) * rSges(9,3) + g(2) * t82) * t70 + (g(1) * (-pkin(15) - t82) + g(2) * rSges(9,3)) * t66) - m(10) * ((g(1) * rSges(10,3) + g(2) * t80) * t70 + (g(2) * rSges(10,3) - g(1) * t80) * t66) - m(11) * ((g(1) * rSges(11,3) - g(2) * t79) * t70 + (g(2) * rSges(11,3) + g(1) * t79) * t66), -m(3) * (g(3) * t91 + t137 * (-rSges(3,1) * t69 + rSges(3,2) * t65)) - m(4) * (g(3) * (t138 - t134) + t137 * (t89 - t133)) - m(5) * (g(1) * (t9 + t77) + g(2) * (t8 + t78) + g(3) * (t96 - t134)) - m(6) * (g(1) * (t9 + t102) + g(2) * (t8 + t103) + g(3) * (t74 - t134) + t72) - m(8) * (g(1) * (t70 * t94 + t30) + g(2) * (t66 * t94 + t29) + g(3) * (-t84 - t134)) - m(11) * (g(1) * (t13 * t70 + t109) + g(2) * (t13 * t66 + t110) + g(3) * (t95 - t134)) + t97, -m(6) * t73 + (-m(4) * t138 - m(5) * t96 - m(6) * t74) * g(3) + t137 * (-m(4) * t89 - m(5) * (t87 - t130) - m(6) * (-t130 + t136)), -m(5) * (g(1) * t77 + g(2) * t78 + g(3) * t139) - m(6) * (g(3) * t76 + t72 + t73), -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t64 - rSges(6,2) * t68) * t47), -m(7) * (g(3) * t86 + t137 * (-rSges(7,1) * t63 - rSges(7,2) * t67)), -m(8) * (g(1) * (-t123 * t70 + t30) + g(2) * (-t123 * t66 + t29) - g(3) * t84) - m(11) * (g(1) * (-t131 * t70 + t109) + g(2) * (-t131 * t66 + t110) + g(3) * t95), t97, -m(10) * (g(1) * (-t119 * t70 + t23) + g(2) * (-t119 * t66 + t20) + g(3) * t108), -m(11) * (g(1) * t109 + g(2) * t110 + g(3) * t98), 0, 0, 0];
taug = t1(:);
