% Calculate potential energy for
% palh1m1DE2
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
% Datum: 2020-04-15 19:16
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE2_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-14 20:02:07
% EndTime: 2020-04-14 20:02:14
% DurationCPUTime: 4.24s
% Computational Cost: add. (70670->211), mult. (105932->266), div. (5052->9), fcn. (67021->48), ass. (0->134)
t167 = pkin(12) + rSges(6,3);
t107 = pkin(1) ^ 2;
t81 = sin(qJ(2));
t83 = sin(pkin(19));
t87 = cos(qJ(2));
t89 = cos(pkin(19));
t57 = t81 * t89 - t83 * t87;
t150 = pkin(7) * t57;
t135 = -0.2e1 * pkin(1) * t150 + t107;
t124 = pkin(7) ^ 2 + t135;
t47 = 0.1e1 / t124;
t138 = 0.1e1 / pkin(3) * t47;
t86 = cos(qJ(3));
t161 = -t86 / 0.2e1;
t156 = -pkin(8) + pkin(3);
t157 = -pkin(8) - pkin(3);
t42 = sqrt(-((pkin(7) - t156) * (pkin(7) + t156) + t135) * ((pkin(7) - t157) * (pkin(7) + t157) + t135));
t58 = t81 * t83 + t87 * t89;
t145 = t42 * t58;
t134 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t45 = t124 - t134;
t50 = pkin(1) - t150;
t38 = -pkin(7) * t145 + t45 * t50;
t165 = -t38 / 0.2e1;
t39 = pkin(7) * t45 * t58 + t42 * t50;
t80 = sin(qJ(3));
t35 = (t39 * t161 + t80 * t165) * t138;
t164 = t39 / 0.2e1;
t36 = (t38 * t161 + t80 * t164) * t138;
t71 = pkin(23) + pkin(22);
t64 = sin(t71);
t65 = cos(t71);
t17 = t35 * t65 + t36 * t64;
t153 = pkin(5) * t17;
t136 = -0.2e1 * pkin(4) * t153 + pkin(5) ^ 2;
t121 = pkin(4) ^ 2 + t136;
t12 = 0.1e1 / t121;
t146 = t12 / pkin(11);
t162 = sin(pkin(21)) / 0.2e1;
t139 = pkin(9) ^ 2 - pkin(11) ^ 2;
t11 = t121 - t139;
t14 = -pkin(4) * t17 + pkin(5);
t18 = -t35 * t64 + t36 * t65;
t154 = pkin(11) - pkin(9);
t155 = -pkin(9) - pkin(11);
t9 = sqrt(-((pkin(4) - t154) * (pkin(4) + t154) + t136) * ((pkin(4) - t155) * (pkin(4) + t155) + t136));
t149 = t18 * t9;
t7 = -pkin(4) * t149 + t11 * t14;
t72 = qJ(2) + qJ(3);
t77 = cos(pkin(21));
t8 = pkin(4) * t11 * t18 + t14 * t9;
t3 = atan2((t7 * t162 + t8 * t77 / 0.2e1) * t146, (-t7 * t77 / 0.2e1 + t8 * t162) * t146) + t72;
t2 = cos(t3);
t166 = t2 * pkin(10);
t163 = sin(pkin(23)) / 0.2e1;
t160 = cos(pkin(18)) / 0.2e1;
t159 = -pkin(2) - pkin(13);
t158 = -pkin(2) + pkin(13);
t75 = sin(pkin(20));
t78 = cos(pkin(20));
t152 = pkin(6) * (-t75 * t86 - t78 * t80);
t151 = pkin(6) * (t75 * t80 - t78 * t86);
t148 = 0.1e1 / pkin(2) / 0.2e1;
t144 = t47 / pkin(8);
t79 = sin(qJ(4));
t82 = sin(qJ(1));
t143 = t79 * t82;
t88 = cos(qJ(1));
t142 = t79 * t88;
t85 = cos(qJ(4));
t141 = t82 * t85;
t140 = t85 * t88;
t76 = cos(pkin(23));
t22 = qJ(2) + atan2((t38 * t163 + t76 * t164) * t138, (t39 * t163 + t76 * t165) * t138);
t100 = pkin(6) ^ 2;
t133 = t100 + t107;
t132 = pkin(1) * t152;
t49 = -0.2e1 * t132;
t122 = t49 + t133;
t120 = 0.1e1 / t122 * t148;
t137 = t100 + t49;
t41 = sqrt(-((pkin(1) - t158) * (pkin(1) + t158) + t137) * ((pkin(1) - t159) * (pkin(1) + t159) + t137));
t105 = pkin(2) ^ 2;
t91 = pkin(13) ^ 2;
t43 = t105 - t91 + t122;
t48 = -pkin(1) + t152;
t33 = qJ(2) + atan2((t43 * t151 - t41 * t48) * t120, (-t41 * t151 - t43 * t48) * t120);
t131 = pkin(14) + r_base(3);
t63 = -t81 * pkin(1) + pkin(16);
t67 = cos(t72);
t59 = pkin(5) * t67 + t63;
t130 = t82 * t59 + r_base(2);
t129 = t88 * t59 + r_base(1);
t128 = t82 * t63 + r_base(2);
t127 = t88 * t63 + r_base(1);
t126 = t82 * pkin(16) + r_base(2);
t125 = t88 * pkin(16) + r_base(1);
t123 = t12 / pkin(9) / 0.2e1;
t119 = 0.1e1 / pkin(13) * t148;
t19 = pkin(22) - t22;
t118 = t87 * pkin(1) + t131;
t66 = sin(t72);
t117 = pkin(5) * t66 + t118;
t1 = sin(t3);
t116 = rSges(5,1) * t2 - rSges(5,2) * t1;
t115 = -rSges(3,1) * t81 - rSges(3,2) * t87;
t114 = rSges(4,1) * t67 - rSges(4,2) * t66;
t20 = sin(t22);
t21 = cos(t22);
t113 = -rSges(8,1) * t20 - rSges(8,2) * t21;
t31 = sin(t33);
t32 = cos(t33);
t112 = -rSges(9,1) * t31 - rSges(9,2) * t32;
t10 = t121 + t139;
t13 = -pkin(4) + t153;
t6 = -atan2((pkin(5) * t10 * t18 - t13 * t9) * t123, (-pkin(5) * t149 - t10 * t13) * t123) + t19;
t4 = sin(t6);
t5 = cos(t6);
t111 = -rSges(11,1) * t4 + rSges(11,2) * t5 + pkin(4) * sin(t19) + t63;
t44 = t124 + t134;
t51 = pkin(1) * t57 - pkin(7);
t37 = -pkin(1) * t145 - t44 * t51;
t40 = pkin(1) * t44 * t58 - t42 * t51;
t84 = sin(pkin(18));
t26 = atan2((t40 * t160 + t37 * t84 / 0.2e1) * t144, (t37 * t160 - t84 * t40 / 0.2e1) * t144);
t23 = sin(t26);
t24 = cos(t26);
t110 = rSges(7,1) * t24 - rSges(7,2) * t23 - pkin(15);
t29 = atan2(t41 * t119, (t105 + t91 + 0.2e1 * t132 - t133) * t119) + t33;
t27 = sin(t29);
t28 = cos(t29);
t109 = -pkin(2) * t31 + rSges(10,1) * t27 + rSges(10,2) * t28 + pkin(16);
t108 = g(1) * r_base(1) + g(2) * r_base(2);
t15 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t88 - rSges(2,2) * t82 + r_base(1)) + g(2) * (rSges(2,1) * t82 + rSges(2,2) * t88 + r_base(2)) + g(3) * (rSges(2,3) + t131)) - m(3) * (g(1) * (rSges(3,3) * t82 + t115 * t88 + t125) + g(2) * (-rSges(3,3) * t88 + t115 * t82 + t126) + g(3) * (rSges(3,1) * t87 - rSges(3,2) * t81 + t131)) - m(4) * (g(1) * (rSges(4,3) * t82 + t114 * t88 + t127) + g(2) * (-rSges(4,3) * t88 + t114 * t82 + t128) + g(3) * (rSges(4,1) * t66 + rSges(4,2) * t67 + t118)) - m(5) * (g(1) * (rSges(5,3) * t82 + t116 * t88 + t129) + g(2) * (-rSges(5,3) * t88 + t116 * t82 + t130) + g(3) * (rSges(5,1) * t1 + rSges(5,2) * t2 + t117)) - m(6) * (g(1) * (t88 * t166 + (t2 * t140 + t143) * rSges(6,1) + (-t2 * t142 + t141) * rSges(6,2) + t129) + g(2) * (t82 * t166 + (t2 * t141 - t142) * rSges(6,1) + (-t2 * t143 - t140) * rSges(6,2) + t130) + g(3) * (-t167 * t2 + t117) + (g(3) * (rSges(6,1) * t85 - rSges(6,2) * t79 + pkin(10)) + (g(1) * t88 + g(2) * t82) * t167) * t1) - m(7) * (g(3) * (rSges(7,1) * t23 + rSges(7,2) * t24 - pkin(17) + t131) + (-g(2) * rSges(7,3) + g(1) * t110) * t88 + (g(1) * rSges(7,3) + g(2) * t110) * t82 + t108) - m(8) * (g(1) * (rSges(8,3) * t82 + t113 * t88 + t127) + g(2) * (-rSges(8,3) * t88 + t113 * t82 + t128) + g(3) * (rSges(8,1) * t21 - rSges(8,2) * t20 + t118)) - m(9) * (g(1) * (rSges(9,3) * t82 + t112 * t88 + t125) + g(2) * (-rSges(9,3) * t88 + t112 * t82 + t126) + g(3) * (rSges(9,1) * t32 - rSges(9,2) * t31 + t131)) - m(10) * (g(3) * (pkin(2) * t32 - rSges(10,1) * t28 + rSges(10,2) * t27 + t131) + (-g(2) * rSges(10,3) + g(1) * t109) * t88 + (g(1) * rSges(10,3) + g(2) * t109) * t82 + t108) - m(11) * (g(3) * (pkin(4) * cos(t19) - t5 * rSges(11,1) - t4 * rSges(11,2) + t118) + (-g(2) * rSges(11,3) + g(1) * t111) * t88 + (g(1) * rSges(11,3) + g(2) * t111) * t82 + t108);
U = t15;
