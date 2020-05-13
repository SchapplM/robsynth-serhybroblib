% Calculate Gravitation load on the joints for
% palh1m2OL
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh1m2OL_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(3,1),zeros(20,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp1: qJ has to be [13x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp1: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m2OL_gravloadJ_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:42
% EndTime: 2020-05-02 21:17:49
% DurationCPUTime: 1.17s
% Computational Cost: add. (594->157), mult. (634->214), div. (0->0), fcn. (428->24), ass. (0->86)
t65 = sin(qJ(5));
t70 = cos(qJ(5));
t28 = rSges(6,1) * t70 - rSges(6,2) * t65;
t135 = pkin(9) + t28;
t133 = pkin(11) + rSges(6,3);
t63 = qJ(2) + qJ(3);
t57 = qJ(4) + t63;
t46 = sin(t57);
t48 = cos(t57);
t134 = t133 * t46 + t135 * t48;
t62 = qJ(2) + qJ(7);
t50 = pkin(19) - t62;
t40 = -qJ(10) + t50;
t33 = sin(t40);
t34 = cos(t40);
t131 = rSges(11,1) * t34 + rSges(11,2) * t33;
t73 = cos(qJ(1));
t130 = g(2) * t73;
t95 = -rSges(11,1) * t33 + rSges(11,2) * t34;
t93 = pkin(4) * sin(t50) + t95;
t129 = t48 * rSges(5,1) - t46 * rSges(5,2);
t68 = sin(qJ(1));
t92 = g(1) * t73 + t68 * g(2);
t66 = sin(qJ(3));
t71 = cos(qJ(3));
t88 = -rSges(4,1) * t66 - rSges(4,2) * t71;
t4 = rSges(3,1) * m(3) + (pkin(1) - t88) * m(4);
t67 = sin(qJ(2));
t72 = cos(qJ(2));
t18 = (rSges(4,1) * t71 - rSges(4,2) * t66) * m(4);
t9 = -rSges(3,2) * m(3) + t18;
t128 = t4 * t67 - t9 * t72;
t41 = pkin(5) * cos(t63);
t125 = pkin(1) * t67;
t44 = pkin(15) - t125;
t12 = t41 + t44;
t126 = t12 * t130;
t124 = pkin(1) * t72;
t123 = pkin(4) * cos(t50);
t122 = pkin(5) * sin(t63);
t121 = pkin(15) * g(1);
t120 = pkin(15) * g(2);
t55 = cos(t62);
t113 = rSges(8,1) * t55;
t52 = sin(t62);
t111 = rSges(8,2) * t52;
t109 = t48 * t68;
t108 = t48 * t73;
t107 = t131 * t68;
t106 = t131 * t73;
t105 = t133 * t109;
t104 = t133 * t108;
t61 = qJ(2) + qJ(8);
t51 = sin(t61);
t54 = cos(t61);
t56 = qJ(9) + t61;
t45 = sin(t56);
t47 = cos(t56);
t99 = -(g(3) * rSges(10,1) - t92 * rSges(10,2)) * t45 + (-t92 * rSges(10,1) - g(3) * rSges(10,2)) * t47;
t100 = -(-(g(3) * rSges(9,1) - t92 * rSges(9,2)) * t51 + (-t92 * rSges(9,1) - g(3) * rSges(9,2)) * t54) * m(9) + m(10) * ((g(3) * t51 + t92 * t54) * pkin(2) + t99);
t98 = -m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3);
t64 = sin(qJ(6));
t69 = cos(qJ(6));
t96 = t69 * rSges(7,1) - rSges(7,2) * t64;
t94 = t41 + t129;
t91 = t68 * g(1) - t130;
t89 = -t113 - t124;
t86 = -rSges(5,1) * t46 - rSges(5,2) * t48;
t27 = -t65 * rSges(6,1) - t70 * rSges(6,2);
t85 = t52 * rSges(8,1) + t55 * rSges(8,2);
t83 = -t44 + t85;
t82 = -t44 - t93;
t81 = t86 * t68;
t80 = t86 * t73;
t77 = -m(2) * rSges(2,1) + t128 + (-m(3) - m(4)) * pkin(15);
t76 = t41 + t134;
t74 = t92 * t135 * t46;
t31 = t73 * t111;
t30 = t68 * t111;
t20 = -pkin(14) + t96;
t19 = -t122 - t124;
t17 = t88 * m(4);
t11 = -t123 - t124;
t7 = t73 * t19;
t6 = t68 * t19;
t1 = [-(t77 * g(1) + t98 * g(2)) * t68 + (-t98 * g(1) + t77 * g(2)) * t73 - m(5) * (t126 + (g(1) * rSges(5,3) + g(2) * t129) * t73 + (g(1) * (-t12 - t129) + g(2) * rSges(5,3)) * t68) - m(6) * (t126 + (-g(1) * t27 + g(2) * t134) * t73 + (g(1) * (-t12 - t134) - g(2) * t27) * t68) - m(7) * (g(1) * (rSges(7,3) * t73 - t20 * t68) + g(2) * (t68 * rSges(7,3) + t20 * t73)) - m(8) * ((g(1) * rSges(8,3) - g(2) * t83) * t73 + (g(2) * rSges(8,3) + g(1) * t83) * t68) - (-(-rSges(9,3) * g(2) + t121) * t68 + (rSges(9,3) * g(1) + t120) * t73 + (rSges(9,1) * t51 + rSges(9,2) * t54) * t91) * m(9) + m(10) * (-(rSges(10,3) * g(2) - t121) * t68 + (-rSges(10,3) * g(1) - t120) * t73 + (-pkin(2) * t51 + t45 * rSges(10,1) + t47 * rSges(10,2)) * t91) - m(11) * ((g(1) * rSges(11,3) - g(2) * t82) * t73 + (g(2) * rSges(11,3) + g(1) * t82) * t68), -m(5) * (g(1) * (t7 + t80) + g(2) * (t6 + t81)) - m(6) * (g(1) * (t7 + t104) + g(2) * (t6 + t105) - t74) - m(8) * (g(1) * (t89 * t73 + t31) + g(2) * (t89 * t68 + t30)) - m(11) * (g(1) * (t11 * t73 + t106) + g(2) * (t11 * t68 + t107)) + t100 + t92 * (t4 * t72 + t67 * t9) + (-m(5) * (t94 - t125) - m(6) * (t76 - t125) - m(8) * (-t85 - t125) - m(11) * (t93 - t125) + t128) * g(3), -m(6) * (g(1) * (-t73 * t122 + t104) + g(2) * (-t68 * t122 + t105) - t74) + (-m(5) * t94 - m(6) * t76 - t17 * t67 - t18 * t72) * g(3) + (-t17 * t72 + t18 * t67 - m(5) * (t86 - t122)) * t92, -m(5) * (g(1) * t80 + g(2) * t81 + g(3) * t129) - m(6) * (g(1) * t104 + g(2) * t105 + g(3) * t134 - t74), -m(6) * (g(1) * (t27 * t108 + t68 * t28) + g(2) * (t27 * t109 - t28 * t73) + g(3) * t27 * t46), -m(7) * (g(3) * t96 + t92 * (-rSges(7,1) * t64 - rSges(7,2) * t69)), -m(8) * (g(1) * (-t73 * t113 + t31) + g(2) * (-t68 * t113 + t30) - g(3) * t85) - m(11) * (g(1) * (-t73 * t123 + t106) + g(2) * (-t68 * t123 + t107) + g(3) * t93), t100, m(10) * t99, -m(11) * (g(1) * t106 + g(2) * t107 + g(3) * t95), 0, 0, 0];
taug = t1(:);
