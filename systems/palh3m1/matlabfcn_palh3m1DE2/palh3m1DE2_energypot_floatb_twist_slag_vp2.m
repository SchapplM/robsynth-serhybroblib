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
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 16:51
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1DE2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1DE2_energypot_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-19 19:34:24
% EndTime: 2020-04-19 19:34:30
% DurationCPUTime: 3.92s
% Computational Cost: add. (69688->132), mult. (104519->152), div. (4968->6), fcn. (66105->38), ass. (0->92)
t62 = sin(qJ(4));
t68 = cos(qJ(4));
t130 = pkin(9) * m(6) + mrSges(6,1) * t68 - mrSges(6,2) * t62 + mrSges(5,1);
t129 = -pkin(11) * m(6) + mrSges(5,2) - mrSges(6,3);
t128 = -m(4) - m(8);
t127 = -m(6) - m(5);
t126 = -m(7) - m(3) - m(2);
t124 = -m(1) - m(9) + t126;
t64 = sin(qJ(2));
t66 = sin(pkin(16));
t70 = cos(qJ(2));
t72 = cos(pkin(16));
t42 = t64 * t66 - t70 * t72;
t110 = pkin(5) * t42;
t103 = (-0.2e1 * t110 + pkin(1)) * pkin(1);
t97 = pkin(5) ^ 2 + t103;
t36 = 0.1e1 / t97;
t105 = t36 / pkin(2);
t117 = sin(qJ(3)) / 0.2e1;
t114 = -pkin(6) + pkin(2);
t115 = -pkin(6) - pkin(2);
t33 = sqrt(-((pkin(5) - t114) * (pkin(5) + t114) + t103) * ((pkin(5) - t115) * (pkin(5) + t115) + t103));
t43 = t64 * t72 + t66 * t70;
t107 = t33 * t43;
t101 = pkin(2) ^ 2 - pkin(6) ^ 2;
t35 = t97 + t101;
t37 = pkin(1) - t110;
t30 = -pkin(5) * t107 + t35 * t37;
t121 = -t30 / 0.2e1;
t31 = pkin(5) * t35 * t43 + t33 * t37;
t69 = cos(qJ(3));
t27 = (t31 * t117 + t121 * t69) * t105;
t120 = t31 / 0.2e1;
t28 = (t117 * t30 + t120 * t69) * t105;
t56 = pkin(18) + pkin(19);
t49 = sin(t56);
t50 = cos(t56);
t18 = -t27 * t50 - t28 * t49;
t111 = pkin(4) * t18;
t104 = -0.2e1 * pkin(3) * t111 + pkin(4) ^ 2;
t98 = pkin(3) ^ 2 + t104;
t12 = 0.1e1 / t98;
t108 = t12 / pkin(10);
t119 = sin(pkin(17)) / 0.2e1;
t57 = qJ(2) + qJ(3);
t59 = cos(pkin(17));
t17 = t27 * t49 - t28 * t50;
t112 = -pkin(8) + pkin(10);
t113 = -pkin(8) - pkin(10);
t9 = sqrt(-((pkin(3) - t112) * (pkin(3) + t112) + t104) * ((pkin(3) - t113) * (pkin(3) + t113) + t104));
t109 = t17 * t9;
t102 = pkin(8) ^ 2 - pkin(10) ^ 2;
t11 = t98 - t102;
t14 = -pkin(3) * t18 + pkin(4);
t7 = -pkin(3) * t109 + t11 * t14;
t8 = pkin(3) * t11 * t17 + t14 * t9;
t3 = atan2((t8 * t59 / 0.2e1 + t7 * t119) * t108, (-t7 * t59 / 0.2e1 + t8 * t119) * t108) + t57;
t1 = sin(t3);
t118 = sin(pkin(19)) / 0.2e1;
t61 = cos(pkin(19));
t22 = qJ(2) + atan2((t118 * t30 + t120 * t61) * t105, (t31 * t118 + t121 * t61) * t105);
t19 = sin(t22);
t2 = cos(t3);
t20 = cos(t22);
t21 = pkin(18) - t22;
t106 = t36 / pkin(6);
t116 = cos(pkin(15)) / 0.2e1;
t34 = t97 - t101;
t38 = pkin(1) * t42 - pkin(5);
t29 = -pkin(1) * t107 - t34 * t38;
t32 = pkin(1) * t34 * t43 - t33 * t38;
t67 = sin(pkin(15));
t26 = atan2((t32 * t116 - t29 * t67 / 0.2e1) * t106, (t29 * t116 + t32 * t67 / 0.2e1) * t106);
t23 = sin(t26);
t24 = cos(t26);
t10 = t98 + t102;
t13 = -pkin(3) + t111;
t96 = t12 / pkin(8) / 0.2e1;
t6 = -atan2((pkin(4) * t10 * t17 - t13 * t9) * t96, (-pkin(4) * t109 - t10 * t13) * t96) + t21;
t4 = sin(t6);
t48 = t70 * pkin(1) + pkin(13);
t5 = cos(t6);
t51 = sin(t57);
t52 = cos(t57);
t123 = m(7) * pkin(7) - m(3) * pkin(13) - t70 * mrSges(3,1) - t24 * mrSges(7,1) + t64 * mrSges(3,2) + t23 * mrSges(7,2) - mrSges(2,1) - m(9) * (pkin(3) * cos(t21) + t48) + t5 * mrSges(9,1) + t4 * mrSges(9,2) - mrSges(8,1) * t20 + mrSges(8,2) * t19 + mrSges(4,1) * t52 - mrSges(4,2) * t51 + t130 * t2 - t129 * t1;
t122 = -mrSges(6,1) * t62 - mrSges(6,2) * t68 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(7,3) - mrSges(8,3) - mrSges(9,3);
t55 = pkin(12) + r_base(3);
t47 = t64 * pkin(1) + t55;
t71 = cos(qJ(1));
t65 = sin(qJ(1));
t44 = -pkin(4) * t52 + t48;
t15 = (m(9) * pkin(3) * sin(t21) - m(1) * r_base(3) - m(7) * pkin(14) - mrSges(3,1) * t64 + mrSges(4,1) * t51 - t23 * mrSges(7,1) - mrSges(8,1) * t19 - t4 * mrSges(9,1) - mrSges(3,2) * t70 + mrSges(4,2) * t52 - t24 * mrSges(7,2) - mrSges(8,2) * t20 + t5 * mrSges(9,2) - mrSges(1,3) - mrSges(2,3) + t127 * (-pkin(4) * t51 + t47) + t126 * t55 + (-m(9) + t128) * t47 + t129 * t2 + t130 * t1) * g(3) + (-mrSges(1,2) + t128 * (t48 * t65 + r_base(2)) + t127 * (t65 * t44 + r_base(2)) + t124 * r_base(2) - t122 * t71 + t123 * t65) * g(2) + (-mrSges(1,1) + t127 * (t71 * t44 + r_base(1)) + t128 * (t48 * t71 + r_base(1)) + t124 * r_base(1) + t123 * t71 + t122 * t65) * g(1);
U = t15;
