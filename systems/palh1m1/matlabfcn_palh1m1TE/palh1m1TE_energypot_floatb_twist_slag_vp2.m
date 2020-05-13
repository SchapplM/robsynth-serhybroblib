% Calculate potential energy for
% palh1m1TE
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
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-13 14:34
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh1m1TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(23,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [23 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp2: pkin has to be [23x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1TE_energypot_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-12 20:10:18
% EndTime: 2020-04-12 20:10:25
% DurationCPUTime: 6.16s
% Computational Cost: add. (80084->192), mult. (120447->262), div. (5724->9), fcn. (76231->28), ass. (0->135)
t88 = pkin(6) ^ 2;
t91 = pkin(1) ^ 2;
t146 = t88 + t91;
t71 = sin(pkin(20));
t73 = cos(pkin(20));
t75 = sin(qJ(3));
t80 = cos(qJ(3));
t156 = pkin(6) * (-t71 * t80 - t73 * t75);
t139 = pkin(1) * t156;
t51 = -0.2e1 * t139;
t90 = 0.1e1 / pkin(2);
t152 = 0.1e1 / (t51 + t146) * t90;
t76 = sin(qJ(2));
t164 = -t76 / 0.2e1;
t155 = pkin(6) * (t71 * t75 - t73 * t80);
t148 = t51 + t88;
t158 = pkin(13) - pkin(2);
t163 = -pkin(2) - pkin(13);
t37 = sqrt(-((pkin(1) - t158) * (pkin(1) + t158) + t148) * ((pkin(1) - t163) * (pkin(1) + t163) + t148));
t134 = -pkin(13) ^ 2 + t146;
t89 = pkin(2) ^ 2;
t40 = t51 + t89 + t134;
t50 = -pkin(1) + t156;
t35 = -t155 * t37 - t40 * t50;
t36 = t155 * t40 - t37 * t50;
t81 = cos(qJ(2));
t120 = (t81 * t35 / 0.2e1 + t36 * t164) * t152;
t149 = 0.1e1 / pkin(13) * t90;
t137 = mrSges(10,2) * t149;
t138 = mrSges(10,1) * t149;
t165 = t134 / 0.2e1 - t89 / 0.2e1 - t139;
t171 = -m(10) * pkin(2) - t137 * t37 / 0.2e1 - t165 * t138 - mrSges(9,1);
t174 = -t165 * t137 + t138 * t37 / 0.2e1 - mrSges(9,2);
t26 = (-t81 * t36 / 0.2e1 + t35 * t164) * t152;
t78 = sin(pkin(19));
t83 = cos(pkin(19));
t59 = t76 * t83 - t78 * t81;
t154 = pkin(7) * t59;
t147 = -0.2e1 * pkin(1) * t154 + t91;
t161 = pkin(7) + pkin(8);
t162 = pkin(7) - pkin(8);
t38 = sqrt(-((-pkin(3) + t161) * (pkin(3) + t162) + t147) * ((pkin(3) + t161) * (-pkin(3) + t162) + t147));
t62 = t76 * t78 + t81 * t83;
t150 = t62 * t38;
t135 = pkin(7) ^ 2 + t147;
t145 = -pkin(3) ^ 2 + pkin(8) ^ 2;
t41 = t135 + t145;
t52 = pkin(1) * t59 - pkin(7);
t125 = -pkin(1) * t150 - t41 * t52;
t126 = pkin(1) * t41 * t62 - t38 * t52;
t133 = cos(pkin(18)) / 0.2e1;
t43 = 0.1e1 / t135;
t151 = t43 / pkin(8);
t153 = sin(pkin(18));
t29 = (t125 * t133 - t153 * t126 / 0.2e1) * t151;
t30 = (t126 * t133 + t125 * t153 / 0.2e1) * t151;
t183 = m(7) * pkin(15) + mrSges(3,1) * t76 - t29 * mrSges(7,1) + mrSges(3,2) * t81 + t30 * mrSges(7,2) - t120 * t174 + t26 * t171 - mrSges(2,1);
t182 = -m(2) - m(7);
t181 = -m(5) - m(6);
t180 = pkin(4) * m(11);
t179 = -m(3) - m(10) - m(9);
t178 = -m(1) + t182;
t177 = -m(11) - m(8) - m(4);
t175 = m(6) * pkin(12) - mrSges(5,2) + mrSges(6,3);
t74 = sin(qJ(4));
t79 = cos(qJ(4));
t173 = -pkin(10) * m(6) - mrSges(6,1) * t79 + mrSges(6,2) * t74 - mrSges(5,1);
t169 = t74 * mrSges(6,1) + t79 * mrSges(6,2) - mrSges(2,2) + mrSges(11,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3) + mrSges(10,3);
t168 = 0.1e1 / pkin(3);
t160 = -pkin(9) - pkin(11);
t159 = pkin(11) - pkin(9);
t157 = pkin(1) * t76;
t144 = -pkin(9) ^ 2 + pkin(11) ^ 2;
t143 = cos(pkin(21));
t142 = cos(pkin(23));
t141 = sin(pkin(21));
t140 = sin(pkin(23));
t69 = pkin(14) + r_base(3);
t77 = sin(qJ(1));
t64 = t77 * pkin(16) + r_base(2);
t82 = cos(qJ(1));
t65 = t82 * pkin(16) + r_base(1);
t136 = pkin(23) + pkin(22);
t132 = pkin(1) - t154;
t63 = t81 * pkin(1) + t69;
t131 = cos(t136);
t130 = sin(t136);
t128 = t135 - t145;
t118 = t168 * (-pkin(7) * t150 + t128 * t132);
t116 = -t118 / 0.2e1;
t119 = t168 * (pkin(7) * t128 * t62 + t132 * t38);
t117 = t119 / 0.2e1;
t27 = (t116 * t142 + t117 * t140) * t43;
t28 = (t142 * t117 + t140 * t118 / 0.2e1) * t43;
t19 = t27 * t81 - t28 * t76;
t20 = -t27 * t76 - t28 * t81;
t60 = t75 * t81 + t76 * t80;
t61 = -t75 * t76 + t80 * t81;
t55 = -t157 * t77 + t64;
t56 = -t157 * t82 + t65;
t115 = t75 * t116 - t80 * t119 / 0.2e1;
t114 = t116 * t80 + t117 * t75;
t113 = t43 * (t114 * t130 + t115 * t131);
t112 = t43 * (t114 * t131 - t115 * t130);
t111 = pkin(5) * t113;
t110 = pkin(4) - t111;
t109 = -pkin(4) * t113 + pkin(5);
t108 = -0.2e1 * pkin(4) * t111 + pkin(5) ^ 2;
t107 = pkin(4) ^ 2 + t108;
t106 = 0.1e1 / t107;
t105 = 0.1e1 / pkin(9) * t106;
t104 = 0.1e1 / pkin(11) * t106;
t103 = t107 + t144;
t102 = t107 - t144;
t101 = sqrt(-((pkin(4) - t159) * (pkin(4) + t159) + t108) * ((pkin(4) - t160) * (pkin(4) + t160) + t108));
t100 = t101 * t112;
t99 = (-pkin(5) * t100 + t102 * t110) * t105;
t98 = (-pkin(4) * t100 + t103 * t109) * t104;
t97 = -(pkin(5) * t102 * t112 + t101 * t110) * t105 / 0.2e1;
t96 = (pkin(4) * t103 * t112 + t101 * t109) * t104 / 0.2e1;
t72 = cos(pkin(22));
t70 = sin(pkin(22));
t49 = t60 * t82;
t48 = t61 * t82;
t47 = t60 * t77;
t46 = t61 * t77;
t18 = t19 * t82;
t17 = t20 * t82;
t16 = t19 * t77;
t15 = t20 * t77;
t10 = -t143 * t98 / 0.2e1 + t141 * t96;
t9 = t141 * t98 / 0.2e1 + t143 * t96;
t8 = -t72 * t99 / 0.2e1 + t70 * t97;
t7 = t70 * t99 / 0.2e1 + t72 * t97;
t1 = (-t30 * mrSges(7,1) - (t19 * t72 - t20 * t70) * t180 - mrSges(8,1) * t19 - mrSges(8,2) * t20 - mrSges(4,1) * t60 - mrSges(4,2) * t61 + m(7) * pkin(17) - mrSges(3,1) * t81 + mrSges(3,2) * t76 - t29 * mrSges(7,2) - m(1) * r_base(3) - (t19 * t8 + t20 * t7) * mrSges(11,1) - (-t19 * t7 + t20 * t8) * mrSges(11,2) - mrSges(1,3) - mrSges(2,3) + t177 * t63 + t173 * (t10 * t60 + t61 * t9) + t181 * (t60 * pkin(5) + t63) + t175 * (t10 * t61 - t60 * t9) + t174 * t26 + t171 * t120 + (t179 + t182) * t69) * g(3) + (-(t15 * t72 + t16 * t70) * t180 - (t15 * t8 - t16 * t7) * mrSges(11,1) - (-t15 * t7 - t16 * t8) * mrSges(11,2) - mrSges(4,1) * t46 + mrSges(4,2) * t47 - mrSges(8,1) * t15 + mrSges(8,2) * t16 - mrSges(1,2) + t179 * t64 + t177 * t55 + t181 * (t46 * pkin(5) + t55) + t175 * (-t10 * t47 - t46 * t9) + t178 * r_base(2) + t173 * (t10 * t46 - t47 * t9)) * g(2) + (-(t17 * t72 + t18 * t70) * t180 - mrSges(4,1) * t48 + mrSges(4,2) * t49 - mrSges(8,1) * t17 + mrSges(8,2) * t18 - (t17 * t8 - t18 * t7) * mrSges(11,1) - (-t17 * t7 - t18 * t8) * mrSges(11,2) - mrSges(1,1) + t179 * t65 + t177 * t56 + t181 * (t48 * pkin(5) + t56) + t175 * (-t10 * t49 - t48 * t9) + t178 * r_base(1) + t173 * (t10 * t48 - t49 * t9)) * g(1) + (g(1) * t183 + t169 * g(2)) * t82 + (-t169 * g(1) + t183 * g(2)) * t77;
U = t1;
