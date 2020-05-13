% Calculate potential energy for
% palh1m1DE1
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
% Datum: 2020-04-14 19:47
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1DE1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp1: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp1: m has to be [11x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [11,3]), ...
  'palh1m1DE1_energypot_floatb_twist_slag_vp1: rSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-13 14:47:36
% EndTime: 2020-04-13 14:47:45
% DurationCPUTime: 7.26s
% Computational Cost: add. (159842->239), mult. (240221->335), div. (11448->9), fcn. (152251->46), ass. (0->144)
t157 = pkin(12) + rSges(6,3);
t125 = pkin(1) ^ 2;
t101 = sin(pkin(19));
t105 = cos(qJ(2));
t107 = cos(pkin(19));
t99 = sin(qJ(2));
t79 = -t101 * t105 + t107 * t99;
t161 = pkin(7) * t79;
t149 = -0.2e1 * pkin(1) * t161 + t125;
t168 = -pkin(8) + pkin(3);
t169 = -pkin(3) - pkin(8);
t59 = sqrt(-((pkin(7) - t168) * (pkin(7) + t168) + t149) * ((pkin(7) - t169) * (pkin(7) + t169) + t149));
t82 = t101 * t99 + t105 * t107;
t155 = t59 * t82;
t141 = pkin(7) ^ 2 + t149;
t146 = pkin(3) ^ 2 - pkin(8) ^ 2;
t62 = t141 + t146;
t73 = pkin(1) - t161;
t52 = -pkin(7) * t155 + t62 * t73;
t175 = -t52 / 0.2e1;
t53 = pkin(7) * t62 * t82 + t59 * t73;
t174 = t53 / 0.2e1;
t173 = sin(pkin(23)) / 0.2e1;
t172 = sin(pkin(21)) / 0.2e1;
t171 = -pkin(2) - pkin(13);
t170 = -pkin(2) + pkin(13);
t167 = -pkin(9) - pkin(11);
t166 = -pkin(9) + pkin(11);
t165 = pkin(1) * t99;
t64 = 0.1e1 / t141;
t152 = 0.1e1 / pkin(3) * t64;
t104 = cos(qJ(3));
t160 = -t104 / 0.2e1;
t98 = sin(qJ(3));
t49 = (t53 * t160 + t98 * t175) * t152;
t50 = (t52 * t160 + t98 * t174) * t152;
t88 = pkin(23) + pkin(22);
t83 = sin(t88);
t84 = cos(t88);
t38 = t49 * t84 + t50 * t83;
t164 = pkin(5) * t38;
t92 = sin(pkin(20));
t96 = cos(pkin(20));
t163 = pkin(6) * (-t104 * t92 - t96 * t98);
t162 = pkin(6) * (-t104 * t96 + t92 * t98);
t159 = cos(pkin(18)) / 0.2e1;
t158 = 0.1e1 / pkin(2) / 0.2e1;
t150 = -0.2e1 * pkin(4) * t164 + pkin(5) ^ 2;
t19 = sqrt(-((pkin(4) - t166) * (pkin(4) + t166) + t150) * ((pkin(4) - t167) * (pkin(4) + t167) + t150));
t39 = -t49 * t83 + t50 * t84;
t156 = t19 * t39;
t140 = pkin(4) ^ 2 + t150;
t34 = 0.1e1 / t140;
t154 = 0.1e1 / pkin(11) * t34;
t153 = 0.1e1 / pkin(8) * t64;
t118 = pkin(6) ^ 2;
t145 = pkin(1) * t163;
t72 = -0.2e1 * t145;
t151 = t118 + t72;
t148 = pkin(9) ^ 2 - pkin(11) ^ 2;
t147 = t118 + t125;
t144 = pkin(14) + r_base(3);
t100 = sin(qJ(1));
t143 = t100 * pkin(16) + r_base(2);
t106 = cos(qJ(1));
t142 = t106 * pkin(16) + r_base(1);
t139 = 0.1e1 / pkin(9) * t34 / 0.2e1;
t138 = 0.1e1 / (t72 + t147) * t158;
t137 = -pkin(13) ^ 2 + t147;
t136 = 0.1e1 / pkin(13) * t158;
t135 = t105 * pkin(1) + t144;
t80 = t104 * t99 + t105 * t98;
t134 = t80 * pkin(5) + t135;
t133 = -rSges(3,1) * t99 - rSges(3,2) * t105;
t93 = cos(pkin(23));
t44 = atan2((t52 * t173 + t93 * t174) * t152, (t53 * t173 + t93 * t175) * t152);
t40 = sin(t44);
t41 = cos(t44);
t25 = t105 * t41 - t40 * t99;
t24 = -t105 * t40 - t41 * t99;
t58 = sqrt(-((pkin(1) - t170) * (pkin(1) + t170) + t151) * ((pkin(1) - t171) * (pkin(1) + t171) + t151));
t123 = pkin(2) ^ 2;
t60 = t123 + t72 + t137;
t71 = -pkin(1) + t163;
t48 = atan2((t60 * t162 - t58 * t71) * t138, (-t58 * t162 - t60 * t71) * t138);
t46 = sin(t48);
t47 = cos(t48);
t33 = t105 * t47 - t46 * t99;
t32 = -t105 * t46 - t47 * t99;
t81 = t104 * t105 - t98 * t99;
t102 = sin(pkin(18));
t61 = t141 - t146;
t74 = pkin(1) * t79 - pkin(7);
t51 = -pkin(1) * t155 - t61 * t74;
t54 = pkin(1) * t61 * t82 - t59 * t74;
t45 = atan2((t54 * t159 + t51 * t102 / 0.2e1) * t153, (t51 * t159 - t102 * t54 / 0.2e1) * t153);
t42 = sin(t45);
t43 = cos(t45);
t132 = rSges(7,1) * t43 - rSges(7,2) * t42 - pkin(15);
t131 = -t100 * t165 + t143;
t130 = -t106 * t165 + t142;
t67 = t81 * t100;
t129 = t67 * pkin(5) + t131;
t69 = t81 * t106;
t128 = t69 * pkin(5) + t130;
t30 = t140 + t148;
t35 = -pkin(4) + t164;
t127 = atan2((pkin(5) * t30 * t39 - t19 * t35) * t139, (-pkin(5) * t156 - t30 * t35) * t139);
t126 = sin(t127);
t103 = cos(qJ(4));
t97 = sin(qJ(4));
t95 = cos(pkin(21));
t94 = cos(pkin(22));
t90 = sin(pkin(22));
t70 = t80 * t106;
t68 = t80 * t100;
t57 = atan2(t58 * t136, (t123 - t137 + 0.2e1 * t145) * t136);
t56 = cos(t57);
t55 = sin(t57);
t36 = -pkin(4) * t38 + pkin(5);
t31 = t140 - t148;
t29 = t32 * t106;
t28 = t33 * t106;
t27 = t32 * t100;
t26 = t33 * t100;
t23 = t24 * t106;
t22 = t25 * t106;
t21 = t24 * t100;
t20 = t25 * t100;
t18 = pkin(4) * t31 * t39 + t19 * t36;
t17 = -pkin(4) * t156 + t31 * t36;
t16 = cos(t127);
t14 = atan2((t17 * t172 + t18 * t95 / 0.2e1) * t154, (-t17 * t95 / 0.2e1 + t18 * t172) * t154);
t13 = cos(t14);
t12 = sin(t14);
t11 = -t90 * t126 - t94 * t16;
t10 = -t94 * t126 + t16 * t90;
t6 = t12 * t81 + t13 * t80;
t5 = t12 * t80 - t81 * t13;
t4 = -t12 * t70 + t13 * t69;
t3 = t12 * t69 + t70 * t13;
t2 = -t12 * t68 + t13 * t67;
t1 = t12 * t67 + t68 * t13;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t106 - rSges(2,2) * t100 + r_base(1)) + g(2) * (rSges(2,1) * t100 + rSges(2,2) * t106 + r_base(2)) + g(3) * (rSges(2,3) + t144)) - m(3) * (g(1) * (rSges(3,3) * t100 + t133 * t106 + t142) + g(2) * (-rSges(3,3) * t106 + t133 * t100 + t143) + g(3) * (rSges(3,1) * t105 - rSges(3,2) * t99 + t144)) - m(4) * (g(1) * (rSges(4,1) * t69 - rSges(4,2) * t70 + rSges(4,3) * t100 + t130) + g(2) * (rSges(4,1) * t67 - rSges(4,2) * t68 - rSges(4,3) * t106 + t131) + g(3) * (rSges(4,1) * t80 + rSges(4,2) * t81 + t135)) - m(5) * (g(1) * (rSges(5,1) * t4 - rSges(5,2) * t3 + rSges(5,3) * t100 + t128) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t1 - rSges(5,3) * t106 + t129) + g(3) * (rSges(5,1) * t6 - rSges(5,2) * t5 + t134)) - m(6) * (g(1) * (t4 * pkin(10) + (t100 * t97 + t103 * t4) * rSges(6,1) + (t100 * t103 - t4 * t97) * rSges(6,2) + t157 * t3 + t128) + g(2) * (t2 * pkin(10) + (t103 * t2 - t106 * t97) * rSges(6,1) + (-t103 * t106 - t2 * t97) * rSges(6,2) + t157 * t1 + t129) + (t134 + (t103 * rSges(6,1) - t97 * rSges(6,2) + pkin(10)) * t6 + t157 * t5) * g(3)) - m(7) * (g(1) * r_base(1) + g(2) * r_base(2) + g(3) * (rSges(7,1) * t42 + rSges(7,2) * t43 - pkin(17) + t144) + (-g(2) * rSges(7,3) + g(1) * t132) * t106 + (g(1) * rSges(7,3) + g(2) * t132) * t100) - m(8) * (g(1) * (rSges(8,1) * t23 - rSges(8,2) * t22 + rSges(8,3) * t100 + t130) + g(2) * (rSges(8,1) * t21 - rSges(8,2) * t20 - rSges(8,3) * t106 + t131) + g(3) * (rSges(8,1) * t25 + rSges(8,2) * t24 + t135)) - m(9) * (g(1) * (rSges(9,1) * t29 - rSges(9,2) * t28 + rSges(9,3) * t100 + t142) + g(2) * (rSges(9,1) * t27 - rSges(9,2) * t26 - rSges(9,3) * t106 + t143) + g(3) * (rSges(9,1) * t33 + rSges(9,2) * t32 + t144)) - m(10) * (g(1) * (t29 * pkin(2) + (t28 * t55 - t29 * t56) * rSges(10,1) + (t28 * t56 + t29 * t55) * rSges(10,2) + t100 * rSges(10,3) + t142) + g(2) * (t27 * pkin(2) + (t26 * t55 - t27 * t56) * rSges(10,1) + (t26 * t56 + t27 * t55) * rSges(10,2) - t106 * rSges(10,3) + t143) + g(3) * (t33 * pkin(2) + (-t32 * t55 - t33 * t56) * rSges(10,1) + (-t32 * t56 + t33 * t55) * rSges(10,2) + t144)) - m(11) * (g(1) * ((-t10 * t22 + t11 * t23) * rSges(11,1) + (-t10 * t23 - t11 * t22) * rSges(11,2) + t100 * rSges(11,3) + t130) + g(2) * ((-t10 * t20 + t11 * t21) * rSges(11,1) + (-t10 * t21 - t11 * t20) * rSges(11,2) - t106 * rSges(11,3) + t131) + g(3) * ((t10 * t24 + t11 * t25) * rSges(11,1) + (-t10 * t25 + t11 * t24) * rSges(11,2) + t135) + (g(1) * (t22 * t90 + t23 * t94) + g(2) * (t20 * t90 + t21 * t94) + g(3) * (-t24 * t90 + t25 * t94)) * pkin(4));
U = t7;
