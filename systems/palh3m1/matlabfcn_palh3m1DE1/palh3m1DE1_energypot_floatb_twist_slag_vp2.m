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
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-19 19:20
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE1_energypot_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-18 10:25:28
% EndTime: 2020-04-18 10:25:39
% DurationCPUTime: 6.72s
% Computational Cost: add. (157276->157), mult. (236358->212), div. (11232->6), fcn. (149713->36), ass. (0->109)
t149 = pkin(3) * m(9);
t148 = -m(2) - m(7);
t147 = -m(6) - m(5);
t146 = -m(1) + t148;
t145 = -m(4) - m(8) - m(9);
t84 = sin(qJ(2));
t86 = sin(pkin(16));
t90 = cos(qJ(2));
t92 = cos(pkin(16));
t63 = t84 * t86 - t90 * t92;
t128 = pkin(5) * t63;
t120 = (-0.2e1 * t128 + pkin(1)) * pkin(1);
t115 = pkin(5) ^ 2 + t120;
t49 = 0.1e1 / t115;
t124 = t49 / pkin(6);
t135 = cos(pkin(15)) / 0.2e1;
t133 = -pkin(6) + pkin(2);
t134 = -pkin(6) - pkin(2);
t46 = sqrt(-((pkin(5) - t133) * (pkin(5) + t133) + t120) * ((pkin(5) - t134) * (pkin(5) + t134) + t120));
t65 = t84 * t92 + t86 * t90;
t125 = t46 * t65;
t119 = -pkin(2) ^ 2 + pkin(6) ^ 2;
t47 = t115 + t119;
t57 = pkin(1) * t63 - pkin(5);
t42 = -pkin(1) * t125 - t47 * t57;
t45 = pkin(1) * t47 * t65 - t46 * t57;
t87 = sin(pkin(15));
t39 = atan2((t45 * t135 - t42 * t87 / 0.2e1) * t124, (t42 * t135 + t45 * t87 / 0.2e1) * t124);
t36 = sin(t39);
t37 = cos(t39);
t144 = m(7) * pkin(7) - t90 * mrSges(3,1) - t37 * mrSges(7,1) + t84 * mrSges(3,2) + t36 * mrSges(7,2) - mrSges(2,1);
t143 = -m(6) * pkin(11) + mrSges(5,2) - mrSges(6,3);
t82 = sin(qJ(4));
t88 = cos(qJ(4));
t142 = -m(6) * pkin(9) - t88 * mrSges(6,1) + t82 * mrSges(6,2) - mrSges(5,1);
t141 = t82 * mrSges(6,1) + t88 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3);
t48 = t115 - t119;
t56 = pkin(1) - t128;
t43 = -pkin(5) * t125 + t48 * t56;
t140 = -t43 / 0.2e1;
t44 = pkin(5) * t48 * t65 + t46 * t56;
t139 = t44 / 0.2e1;
t138 = sin(pkin(17)) / 0.2e1;
t137 = sin(pkin(19)) / 0.2e1;
t83 = sin(qJ(3));
t136 = t83 / 0.2e1;
t132 = -pkin(8) - pkin(10);
t131 = -pkin(8) + pkin(10);
t130 = pkin(1) * t90;
t122 = 0.1e1 / pkin(2) * t49;
t89 = cos(qJ(3));
t40 = (t44 * t136 + t140 * t89) * t122;
t41 = (t136 * t43 + t139 * t89) * t122;
t75 = pkin(18) + pkin(19);
t69 = sin(t75);
t70 = cos(t75);
t33 = -t40 * t70 - t41 * t69;
t129 = pkin(4) * t33;
t121 = -0.2e1 * pkin(3) * t129 + pkin(4) ^ 2;
t19 = sqrt(-((pkin(3) - t131) * (pkin(3) + t131) + t121) * ((pkin(3) - t132) * (pkin(3) + t132) + t121));
t32 = t40 * t69 - t41 * t70;
t127 = t19 * t32;
t114 = pkin(3) ^ 2 + t121;
t28 = 0.1e1 / t114;
t126 = t28 / pkin(10);
t123 = pkin(8) ^ 2 - pkin(10) ^ 2;
t74 = pkin(12) + r_base(3);
t85 = sin(qJ(1));
t118 = t85 * pkin(13) + r_base(2);
t91 = cos(qJ(1));
t117 = t91 * pkin(13) + r_base(1);
t116 = t28 / pkin(8) / 0.2e1;
t60 = t85 * t130 + t118;
t61 = t91 * t130 + t117;
t66 = t84 * pkin(1) + t74;
t80 = cos(pkin(19));
t38 = atan2((t137 * t43 + t139 * t80) * t122, (t44 * t137 + t140 * t80) * t122);
t34 = sin(t38);
t35 = cos(t38);
t25 = t34 * t90 + t35 * t84;
t24 = -t34 * t84 + t35 * t90;
t109 = t83 * t90 + t84 * t89;
t62 = t83 * t84 - t89 * t90;
t26 = t114 + t123;
t29 = -pkin(3) + t129;
t107 = atan2((pkin(4) * t26 * t32 - t19 * t29) * t116, (-pkin(4) * t127 - t26 * t29) * t116);
t106 = sin(t107);
t81 = cos(pkin(18));
t79 = sin(pkin(18));
t77 = cos(pkin(17));
t55 = t62 * t91;
t54 = t109 * t91;
t53 = t62 * t85;
t52 = t109 * t85;
t30 = -pkin(3) * t33 + pkin(4);
t27 = t114 - t123;
t23 = t24 * t91;
t22 = t25 * t91;
t21 = t24 * t85;
t20 = t25 * t85;
t18 = pkin(3) * t27 * t32 + t19 * t30;
t17 = -pkin(3) * t127 + t27 * t30;
t16 = cos(t107);
t14 = atan2((t18 * t77 / 0.2e1 + t17 * t138) * t126, (-t17 * t77 / 0.2e1 + t18 * t138) * t126);
t13 = cos(t14);
t12 = sin(t14);
t11 = -t106 * t79 - t81 * t16;
t10 = -t106 * t81 + t16 * t79;
t1 = (-t37 * mrSges(7,2) - t62 * mrSges(4,2) + t109 * mrSges(4,1) - m(7) * pkin(14) - m(1) * r_base(3) - t36 * mrSges(7,1) - (-t24 * t79 + t25 * t81) * t149 - t24 * mrSges(8,2) - (t10 * t24 + t11 * t25) * mrSges(9,1) - (-t10 * t25 + t11 * t24) * mrSges(9,2) - t25 * mrSges(8,1) - t84 * mrSges(3,1) - t90 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(3) + t148) * t74 + t145 * t66 + t142 * (-t109 * t13 + t12 * t62) + t147 * (-pkin(4) * t109 + t66) + t143 * (-t109 * t12 - t62 * t13)) * g(3) + (-t52 * mrSges(4,2) - t53 * mrSges(4,1) - m(3) * t118 + t20 * mrSges(8,2) - (-t10 * t20 + t11 * t21) * mrSges(9,1) - (-t10 * t21 - t11 * t20) * mrSges(9,2) - t21 * mrSges(8,1) - (t20 * t79 + t21 * t81) * t149 - mrSges(1,2) + t145 * t60 + t147 * (t53 * pkin(4) + t60) + t143 * (t12 * t53 - t52 * t13) + t146 * r_base(2) + t144 * t85 + t142 * (t12 * t52 + t13 * t53) + t141 * t91) * g(2) + (-t54 * mrSges(4,2) - t55 * mrSges(4,1) - m(3) * t117 - (t22 * t79 + t23 * t81) * t149 + t22 * mrSges(8,2) - (-t10 * t22 + t11 * t23) * mrSges(9,1) - (-t10 * t23 - t11 * t22) * mrSges(9,2) - t23 * mrSges(8,1) - mrSges(1,1) + t145 * t61 + t147 * (t55 * pkin(4) + t61) + t143 * (t12 * t55 - t54 * t13) + t146 * r_base(1) + t144 * t91 + t142 * (t12 * t54 + t13 * t55) - t141 * t85) * g(1);
U = t1;
