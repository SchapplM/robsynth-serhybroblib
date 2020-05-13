% Calculate potential energy for
% palh3m1DE1
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
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:25:32
% EndTime: 2020-04-18 10:25:42
% DurationCPUTime: 6.62s
% Computational Cost: add. (157276->185), mult. (236340->263), div. (11232->6), fcn. (149713->36), ass. (0->112)
t128 = pkin(11) + rSges(6,3);
t80 = sin(qJ(2));
t82 = sin(pkin(16));
t86 = cos(qJ(2));
t88 = cos(pkin(16));
t61 = t80 * t82 - t86 * t88;
t129 = pkin(5) * t61;
t119 = (-0.2e1 * t129 + pkin(1)) * pkin(1);
t134 = -pkin(6) + pkin(2);
t135 = -pkin(6) - pkin(2);
t46 = sqrt(-((pkin(5) - t134) * (pkin(5) + t134) + t119) * ((pkin(5) - t135) * (pkin(5) + t135) + t119));
t63 = t80 * t88 + t82 * t86;
t125 = t46 * t63;
t114 = pkin(5) ^ 2 + t119;
t121 = pkin(2) ^ 2 - pkin(6) ^ 2;
t48 = t114 + t121;
t56 = pkin(1) - t129;
t43 = -pkin(5) * t125 + t48 * t56;
t141 = -t43 / 0.2e1;
t44 = pkin(5) * t48 * t63 + t46 * t56;
t140 = t44 / 0.2e1;
t139 = sin(pkin(17)) / 0.2e1;
t138 = sin(pkin(19)) / 0.2e1;
t79 = sin(qJ(3));
t137 = t79 / 0.2e1;
t136 = cos(pkin(15)) / 0.2e1;
t133 = -pkin(8) - pkin(10);
t132 = -pkin(8) + pkin(10);
t131 = pkin(1) * t86;
t49 = 0.1e1 / t114;
t120 = 0.1e1 / pkin(2) * t49;
t85 = cos(qJ(3));
t40 = (t44 * t137 + t141 * t85) * t120;
t41 = (t137 * t43 + t140 * t85) * t120;
t71 = pkin(18) + pkin(19);
t66 = sin(t71);
t67 = cos(t71);
t33 = -t40 * t67 - t41 * t66;
t130 = pkin(4) * t33;
t123 = -0.2e1 * pkin(3) * t130 + pkin(4) ^ 2;
t19 = sqrt(-((pkin(3) - t132) * (pkin(3) + t132) + t123) * ((pkin(3) - t133) * (pkin(3) + t133) + t123));
t32 = t40 * t66 - t41 * t67;
t127 = t19 * t32;
t115 = pkin(3) ^ 2 + t123;
t28 = 0.1e1 / t115;
t126 = t28 / pkin(10);
t124 = t49 / pkin(6);
t122 = pkin(8) ^ 2 - pkin(10) ^ 2;
t118 = pkin(12) + r_base(3);
t81 = sin(qJ(1));
t117 = t81 * pkin(13) + r_base(2);
t87 = cos(qJ(1));
t116 = t87 * pkin(13) + r_base(1);
t113 = t28 / pkin(8) / 0.2e1;
t112 = t80 * pkin(1) + t118;
t111 = t81 * t131 + t117;
t110 = t87 * t131 + t116;
t105 = t79 * t86 + t80 * t85;
t109 = -pkin(4) * t105 + t112;
t60 = t79 * t80 - t85 * t86;
t53 = t60 * t81;
t108 = t53 * pkin(4) + t111;
t55 = t60 * t87;
t107 = t55 * pkin(4) + t110;
t106 = rSges(3,1) * t86 - rSges(3,2) * t80;
t76 = cos(pkin(19));
t38 = atan2((t138 * t43 + t140 * t76) * t120, (t44 * t138 + t141 * t76) * t120);
t34 = sin(t38);
t35 = cos(t38);
t25 = t34 * t86 + t35 * t80;
t24 = -t34 * t80 + t35 * t86;
t47 = t114 - t121;
t57 = pkin(1) * t61 - pkin(5);
t42 = -pkin(1) * t125 - t47 * t57;
t45 = pkin(1) * t47 * t63 - t46 * t57;
t83 = sin(pkin(15));
t39 = atan2((t45 * t136 - t42 * t83 / 0.2e1) * t124, (t42 * t136 + t45 * t83 / 0.2e1) * t124);
t36 = sin(t39);
t37 = cos(t39);
t104 = rSges(7,1) * t37 - rSges(7,2) * t36 - pkin(7);
t26 = t115 + t122;
t29 = -pkin(3) + t130;
t103 = atan2((pkin(4) * t26 * t32 - t19 * t29) * t113, (-pkin(4) * t127 - t26 * t29) * t113);
t102 = sin(t103);
t84 = cos(qJ(4));
t78 = sin(qJ(4));
t77 = cos(pkin(18));
t75 = sin(pkin(18));
t73 = cos(pkin(17));
t54 = t105 * t87;
t52 = t105 * t81;
t30 = -pkin(3) * t33 + pkin(4);
t27 = t115 - t122;
t23 = t24 * t87;
t22 = t25 * t87;
t21 = t24 * t81;
t20 = t25 * t81;
t18 = pkin(3) * t27 * t32 + t19 * t30;
t17 = -pkin(3) * t127 + t27 * t30;
t16 = cos(t103);
t14 = atan2((t18 * t73 / 0.2e1 + t17 * t139) * t126, (-t17 * t73 / 0.2e1 + t18 * t139) * t126);
t13 = cos(t14);
t12 = sin(t14);
t11 = -t102 * t75 - t77 * t16;
t10 = -t102 * t77 + t16 * t75;
t6 = -t105 * t13 + t12 * t60;
t5 = -t105 * t12 - t60 * t13;
t4 = t12 * t54 + t13 * t55;
t3 = t12 * t55 - t54 * t13;
t2 = t12 * t52 + t13 * t53;
t1 = t12 * t53 - t52 * t13;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t87 - rSges(2,2) * t81 + r_base(1)) + g(2) * (rSges(2,1) * t81 + rSges(2,2) * t87 + r_base(2)) + g(3) * (rSges(2,3) + t118)) - m(3) * (g(1) * (rSges(3,3) * t81 + t106 * t87 + t116) + g(2) * (-rSges(3,3) * t87 + t106 * t81 + t117) + g(3) * (rSges(3,1) * t80 + rSges(3,2) * t86 + t118)) - m(4) * (g(1) * (rSges(4,1) * t55 + rSges(4,2) * t54 + rSges(4,3) * t81 + t110) + g(2) * (rSges(4,1) * t53 + rSges(4,2) * t52 - rSges(4,3) * t87 + t111) + g(3) * (-rSges(4,1) * t105 + rSges(4,2) * t60 + t112)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t81 + t107) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t1 - rSges(5,3) * t87 + t108) + g(3) * (rSges(5,1) * t6 - rSges(5,2) * t5 + t109)) - m(6) * (g(1) * (t4 * pkin(9) + (t4 * t84 + t78 * t81) * rSges(6,1) + (-t4 * t78 + t81 * t84) * rSges(6,2) + t128 * t3 + t107) + g(2) * (t2 * pkin(9) + (t2 * t84 - t78 * t87) * rSges(6,1) + (-t2 * t78 - t84 * t87) * rSges(6,2) + t128 * t1 + t108) + (t109 + (t84 * rSges(6,1) - t78 * rSges(6,2) + pkin(9)) * t6 + t128 * t5) * g(3)) - m(7) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(7,1) * t36 + rSges(7,2) * t37 + pkin(14) + t118) + (-g(2) * rSges(7,3) + g(1) * t104) * t87 + (g(1) * rSges(7,3) + g(2) * t104) * t81) - m(8) * (g(1) * (rSges(8,1) * t23 - rSges(8,2) * t22 + rSges(8,3) * t81 + t110) + g(2) * (rSges(8,1) * t21 - rSges(8,2) * t20 - rSges(8,3) * t87 + t111) + g(3) * (rSges(8,1) * t25 + rSges(8,2) * t24 + t112)) - m(9) * (g(1) * ((-t10 * t22 + t11 * t23) * rSges(9,1) + (-t10 * t23 - t11 * t22) * rSges(9,2) + t81 * rSges(9,3) + t110) + g(2) * ((-t10 * t20 + t11 * t21) * rSges(9,1) + (-t10 * t21 - t11 * t20) * rSges(9,2) - t87 * rSges(9,3) + t111) + g(3) * ((t10 * t24 + t11 * t25) * rSges(9,1) + (-t10 * t25 + t11 * t24) * rSges(9,2) + t112) + (g(1) * (t22 * t75 + t23 * t77) + g(2) * (t20 * t75 + t21 * t77) + g(3) * (-t24 * t75 + t25 * t77)) * pkin(3));
U = t7;
