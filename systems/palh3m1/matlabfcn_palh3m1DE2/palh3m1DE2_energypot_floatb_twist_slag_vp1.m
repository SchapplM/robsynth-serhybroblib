% Calculate potential energy for
% palh3m1DE2
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
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp1: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp1: rSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:34:27
% EndTime: 2020-04-19 19:34:33
% DurationCPUTime: 3.95s
% Computational Cost: add. (69688->164), mult. (104501->213), div. (4968->6), fcn. (66105->38), ass. (0->102)
t128 = pkin(11) + rSges(6,3);
t62 = sin(qJ(2));
t64 = sin(pkin(16));
t68 = cos(qJ(2));
t70 = cos(pkin(16));
t42 = t62 * t64 - t68 * t70;
t115 = pkin(5) * t42;
t103 = (-0.2e1 * t115 + pkin(1)) * pkin(1);
t94 = pkin(5) ^ 2 + t103;
t36 = 0.1e1 / t94;
t109 = t36 / pkin(2);
t122 = sin(qJ(3)) / 0.2e1;
t119 = -pkin(6) + pkin(2);
t120 = -pkin(6) - pkin(2);
t33 = sqrt(-((pkin(5) - t119) * (pkin(5) + t119) + t103) * ((pkin(5) - t120) * (pkin(5) + t120) + t103));
t43 = t62 * t70 + t64 * t68;
t111 = t33 * t43;
t101 = pkin(2) ^ 2 - pkin(6) ^ 2;
t35 = t94 + t101;
t37 = pkin(1) - t115;
t30 = -pkin(5) * t111 + t35 * t37;
t126 = -t30 / 0.2e1;
t31 = pkin(5) * t35 * t43 + t33 * t37;
t67 = cos(qJ(3));
t27 = (t31 * t122 + t67 * t126) * t109;
t125 = t31 / 0.2e1;
t28 = (t30 * t122 + t67 * t125) * t109;
t54 = pkin(18) + pkin(19);
t48 = sin(t54);
t49 = cos(t54);
t18 = -t27 * t49 - t28 * t48;
t116 = pkin(4) * t18;
t104 = -0.2e1 * pkin(3) * t116 + pkin(4) ^ 2;
t95 = pkin(3) ^ 2 + t104;
t12 = 0.1e1 / t95;
t112 = t12 / pkin(10);
t124 = sin(pkin(17)) / 0.2e1;
t55 = qJ(2) + qJ(3);
t57 = cos(pkin(17));
t102 = pkin(8) ^ 2 - pkin(10) ^ 2;
t11 = t95 - t102;
t17 = t27 * t48 - t28 * t49;
t117 = -pkin(8) + pkin(10);
t118 = -pkin(8) - pkin(10);
t9 = sqrt(-((pkin(3) - t117) * (pkin(3) + t117) + t104) * ((pkin(3) - t118) * (pkin(3) + t118) + t104));
t114 = t17 * t9;
t14 = -pkin(3) * t18 + pkin(4);
t7 = -pkin(3) * t114 + t11 * t14;
t8 = pkin(3) * t11 * t17 + t14 * t9;
t3 = atan2((t8 * t57 / 0.2e1 + t7 * t124) * t112, (-t7 * t57 / 0.2e1 + t8 * t124) * t112) + t55;
t2 = cos(t3);
t127 = t2 * pkin(9);
t123 = sin(pkin(19)) / 0.2e1;
t121 = cos(pkin(15)) / 0.2e1;
t110 = t36 / pkin(6);
t60 = sin(qJ(4));
t63 = sin(qJ(1));
t108 = t60 * t63;
t69 = cos(qJ(1));
t107 = t60 * t69;
t66 = cos(qJ(4));
t106 = t63 * t66;
t105 = t66 * t69;
t47 = t68 * pkin(1) + pkin(13);
t59 = cos(pkin(19));
t22 = qJ(2) + atan2((t30 * t123 + t59 * t125) * t109, (t31 * t123 + t59 * t126) * t109);
t100 = pkin(12) + r_base(3);
t51 = cos(t55);
t44 = -pkin(4) * t51 + t47;
t99 = t63 * t44 + r_base(2);
t98 = t69 * t44 + r_base(1);
t97 = t63 * t47 + r_base(2);
t96 = t69 * t47 + r_base(1);
t93 = t12 / pkin(8) / 0.2e1;
t21 = pkin(18) - t22;
t92 = t62 * pkin(1) + t100;
t1 = sin(t3);
t91 = -rSges(5,1) * t2 + rSges(5,2) * t1;
t50 = sin(t55);
t90 = -rSges(4,1) * t51 + rSges(4,2) * t50;
t19 = sin(t22);
t20 = cos(t22);
t89 = rSges(8,1) * t20 - rSges(8,2) * t19;
t10 = t95 + t102;
t13 = -pkin(3) + t116;
t6 = -atan2((pkin(4) * t10 * t17 - t13 * t9) * t93, (-pkin(4) * t114 - t10 * t13) * t93) + t21;
t4 = sin(t6);
t5 = cos(t6);
t88 = -rSges(9,1) * t5 - rSges(9,2) * t4 + pkin(3) * cos(t21) + t47;
t34 = t94 - t101;
t38 = pkin(1) * t42 - pkin(5);
t29 = -pkin(1) * t111 - t34 * t38;
t32 = pkin(1) * t34 * t43 - t33 * t38;
t65 = sin(pkin(15));
t26 = atan2((t32 * t121 - t29 * t65 / 0.2e1) * t110, (t29 * t121 + t32 * t65 / 0.2e1) * t110);
t23 = sin(t26);
t24 = cos(t26);
t87 = rSges(7,1) * t24 - rSges(7,2) * t23 - pkin(7);
t86 = rSges(3,1) * t68 - rSges(3,2) * t62 + pkin(13);
t85 = -pkin(4) * t50 + t92;
t84 = g(1) * r_base(1) + g(2) * r_base(2);
t15 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t69 - rSges(2,2) * t63 + r_base(1)) + g(2) * (rSges(2,1) * t63 + rSges(2,2) * t69 + r_base(2)) + g(3) * (rSges(2,3) + t100)) - m(3) * (g(3) * (rSges(3,1) * t62 + rSges(3,2) * t68 + t100) + (-g(2) * rSges(3,3) + g(1) * t86) * t69 + (g(1) * rSges(3,3) + g(2) * t86) * t63 + t84) - m(4) * (g(1) * (rSges(4,3) * t63 + t90 * t69 + t96) + g(2) * (-rSges(4,3) * t69 + t90 * t63 + t97) + g(3) * (-rSges(4,1) * t50 - rSges(4,2) * t51 + t92)) - m(5) * (g(1) * (rSges(5,3) * t63 + t91 * t69 + t98) + g(2) * (-rSges(5,3) * t69 + t91 * t63 + t99) + g(3) * (-rSges(5,1) * t1 - rSges(5,2) * t2 + t85)) - m(6) * (g(1) * (-t69 * t127 + (-t2 * t105 + t108) * rSges(6,1) + (t2 * t107 + t106) * rSges(6,2) + t98) + g(2) * (-t63 * t127 + (-t2 * t106 - t107) * rSges(6,1) + (t2 * t108 - t105) * rSges(6,2) + t99) + g(3) * (t128 * t2 + t85) + (g(3) * (-rSges(6,1) * t66 + rSges(6,2) * t60 - pkin(9)) - (g(1) * t69 + g(2) * t63) * t128) * t1) - m(7) * (g(3) * (rSges(7,1) * t23 + rSges(7,2) * t24 + pkin(14) + t100) + (-g(2) * rSges(7,3) + g(1) * t87) * t69 + (g(1) * rSges(7,3) + g(2) * t87) * t63 + t84) - m(8) * (g(1) * (rSges(8,3) * t63 + t89 * t69 + t96) + g(2) * (-rSges(8,3) * t69 + t89 * t63 + t97) + g(3) * (rSges(8,1) * t19 + rSges(8,2) * t20 + t92)) - m(9) * (g(3) * (-pkin(3) * sin(t21) + t4 * rSges(9,1) - t5 * rSges(9,2) + t92) + (-g(2) * rSges(9,3) + g(1) * t88) * t69 + (g(1) * rSges(9,3) + g(2) * t88) * t63 + t84);
U = t15;
