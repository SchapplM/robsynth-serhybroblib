% Calculate Gravitation load on the joints for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh3m2DE2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(18,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2DE2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:21:08
% EndTime: 2020-05-07 04:21:08
% DurationCPUTime: 0.68s
% Computational Cost: add. (1059->142), mult. (474->150), div. (0->0), fcn. (286->63), ass. (0->88)
t100 = mrSges(9,2) * g(2);
t113 = pkin(15) + qJ(2);
t90 = pkin(18) + t113;
t52 = pkin(17) + qJ(3) + t90;
t31 = atan2(-sin(t52), -cos(t52));
t41 = atan2(sin(t90), -cos(t90));
t108 = pkin(17) - qJ(2) - t41;
t32 = qJ(1) + t108;
t11 = -t31 + t32;
t121 = mrSges(9,2) * g(1);
t129 = mrSges(4,2) * g(1);
t131 = mrSges(9,1) * g(2);
t132 = mrSges(9,1) * g(1);
t103 = m(5) + m(6);
t145 = pkin(4) * t103 + mrSges(4,1);
t152 = g(2) * t145;
t95 = qJ(2) + qJ(3);
t93 = qJ(1) - t95;
t156 = -(t152 + t129) * cos(t93) / 0.2e1 - (t121 + t131) * cos(t11) / 0.2e1 + (-t100 + t132) * sin(t11) / 0.2e1;
t33 = -qJ(1) + t108;
t12 = -t31 + t33;
t92 = qJ(1) + t95;
t154 = (t100 + t132) * sin(t12) / 0.2e1 + (-t121 + t131) * cos(t12) / 0.2e1 + (t152 - t129) * cos(t92) / 0.2e1;
t153 = g(1) * t145;
t151 = t145 * cos(t95) - sin(t95) * mrSges(4,2);
t126 = mrSges(6,2) * g(2);
t127 = mrSges(6,2) * g(1);
t136 = mrSges(6,1) * g(2);
t137 = mrSges(6,1) * g(1);
t47 = pkin(16) + t52;
t30 = atan2(-sin(t47), cos(t47));
t109 = t30 + t95;
t27 = -qJ(4) + t109;
t24 = qJ(1) + t27;
t26 = qJ(4) + t109;
t25 = qJ(1) - t26;
t97 = qJ(1) - qJ(4);
t150 = (cos(t97) / 0.2e1 - cos(t24) / 0.4e1 - cos(t25) / 0.4e1) * (t127 + t136) + (-sin(t97) / 0.2e1 + sin(t25) / 0.4e1 + sin(t24) / 0.4e1) * (-t126 + t137);
t144 = 2 * g(1);
t143 = 2 * pkin(8) * m(6) + 2 * mrSges(5,1);
t29 = -t30 + t93;
t138 = cos(t29);
t135 = mrSges(7,1) * g(1);
t134 = mrSges(7,1) * g(2);
t133 = mrSges(8,1) * g(2);
t130 = mrSges(3,2) * g(1);
t128 = mrSges(4,2) * g(2);
t125 = mrSges(7,2) * g(1);
t124 = mrSges(7,2) * g(2);
t123 = mrSges(8,2) * g(1);
t122 = mrSges(8,2) * g(2);
t80 = m(8) + m(9) + m(4) + t103;
t48 = (pkin(1) * t80 + mrSges(3,1));
t120 = g(2) * t48;
t79 = pkin(10) * m(6) - mrSges(5,2) + mrSges(6,3);
t118 = g(2) * t79;
t117 = t48 * g(1);
t99 = qJ(1) - qJ(2);
t98 = qJ(1) + qJ(2);
t28 = t30 + t92;
t112 = cos(t28) / 0.4e1;
t91 = pkin(14) - t113;
t22 = qJ(1) + t26;
t23 = qJ(1) - t27;
t67 = t126 + t137;
t70 = -t127 + t136;
t96 = qJ(1) + qJ(4);
t107 = t70 * cos(t96) / 0.2e1 - sin(t96) * t67 / 0.2e1 - (sin(t22) + sin(t23)) * t67 / 0.4e1 + (cos(t22) + cos(t23)) * t70 / 0.4e1;
t13 = -t31 + t108;
t71 = sin(t92);
t72 = sin(t93);
t106 = (-t128 + t153) * t72 / 0.2e1 + (-t128 - t153) * t71 / 0.2e1 + (mrSges(9,1) * cos(t13) + mrSges(9,2) * sin(t13)) * g(3) + t154 + t156;
t77 = qJ(1) + t91;
t105 = (t120 + t130) * cos(t99) / 0.2e1 - (-t124 + t135) * sin(t77) / 0.2e1 + (t125 + t134) * cos(t77) / 0.2e1;
t78 = -qJ(1) + t91;
t104 = (-t120 + t130) * cos(t98) / 0.2e1 + (-t124 - t135) * sin(t78) / 0.2e1 + (t125 - t134) * cos(t78) / 0.2e1;
t102 = mrSges(8,1) * g(1);
t101 = mrSges(3,2) * g(2);
t85 = sin(t99);
t84 = sin(t98);
t51 = g(1) * t143;
t50 = (mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3));
t46 = t145 * t144;
t45 = 2 * t117;
t42 = (m(7) * pkin(6) - mrSges(2,1) + (-m(3) - t80) * pkin(12));
t39 = -t41 + t99;
t38 = t41 + t98;
t1 = [t154 - t150 + (g(1) * t50 + g(2) * t42) * cos(qJ(1)) + (-t46 + 0.2e1 * t128) * t72 / 0.4e1 + (-t46 - 0.2e1 * t128) * t71 / 0.4e1 + (-t123 - t133) * cos(t39) / 0.2e1 + (t123 - t133) * cos(t38) / 0.2e1 + (t102 - t122) * sin(t39) / 0.2e1 + (t102 + t122) * sin(t38) / 0.2e1 + (-t51 + 2 * t118) * sin(t28) / 0.4e1 + (-t51 - 2 * t118) * sin(t29) / 0.4e1 + t107 + t104 - t105 + ((-cos(t32) / 0.2e1 - cos(t33) / 0.2e1) * g(2) + (sin(t32) / 0.2e1 - sin(t33) / 0.2e1) * g(1)) * m(9) * pkin(3) - t156 + (t45 + 2 * t101) * t84 / 0.4e1 + (t45 - 2 * t101) * t85 / 0.4e1 - (g(1) * t42 - g(2) * t50) * sin(qJ(1)) + (t112 + t138 / 0.4e1) * g(2) * t143 + (t112 - t138 / 0.4e1) * t79 * t144, t106 + t104 + t105 + (-mrSges(7,1) * cos(t91) + sin(qJ(2)) * mrSges(3,2) - mrSges(7,2) * sin(t91) - t48 * cos(qJ(2)) + t151) * g(3) + (t101 - t117) * t85 / 0.2e1 + (t101 + t117) * t84 / 0.2e1, g(3) * t151 + t106, ((-sin(t27) / 0.2e1 - sin(t26) / 0.2e1) * mrSges(6,2) + (-cos(t27) / 0.2e1 + cos(t26) / 0.2e1) * mrSges(6,1)) * g(3) + t107 + t150];
taug = t1(:);
