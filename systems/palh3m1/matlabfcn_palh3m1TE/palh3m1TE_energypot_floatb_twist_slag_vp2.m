% Calculate potential energy for
% palh3m1TE
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
% Datum: 2020-04-18 10:11
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh3m1TE_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(19,1),zeros(9,1),zeros(9,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [19 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp2: pkin has to be [19x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1TE_energypot_floatb_twist_slag_vp2: mrSges has to be [9x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-17 15:18:46
% EndTime: 2020-04-17 15:18:54
% DurationCPUTime: 5.02s
% Computational Cost: add. (78772->157), mult. (118446->215), div. (5616->6), fcn. (74977->24), ass. (0->104)
t139 = pkin(3) * m(9);
t138 = -m(2) - m(7);
t137 = -m(5) - m(6);
t136 = -m(1) + t138;
t135 = -m(4) - m(9) - m(8);
t109 = cos(pkin(15)) / 0.2e1;
t58 = sin(qJ(2));
t60 = sin(pkin(16));
t63 = cos(qJ(2));
t65 = cos(pkin(16));
t44 = t58 * t60 - t63 * t65;
t124 = pkin(5) * t44;
t120 = (-0.2e1 * t124 + pkin(1)) * pkin(1);
t110 = pkin(5) ^ 2 + t120;
t31 = 0.1e1 / t110;
t122 = t31 / pkin(6);
t123 = sin(pkin(15));
t128 = -pkin(6) + pkin(2);
t129 = -pkin(6) - pkin(2);
t29 = sqrt(-((pkin(5) - t128) * (pkin(5) + t128) + t120) * ((pkin(5) - t129) * (pkin(5) + t129) + t120));
t46 = t58 * t65 + t60 * t63;
t121 = t46 * t29;
t119 = -pkin(2) ^ 2 + pkin(6) ^ 2;
t30 = t110 + t119;
t38 = pkin(1) * t44 - pkin(5);
t98 = -pkin(1) * t121 - t30 * t38;
t99 = pkin(1) * t30 * t46 - t29 * t38;
t23 = (t98 * t109 + t99 * t123 / 0.2e1) * t122;
t24 = (t99 * t109 - t98 * t123 / 0.2e1) * t122;
t134 = m(7) * pkin(7) - t63 * mrSges(3,1) - t23 * mrSges(7,1) + t58 * mrSges(3,2) + t24 * mrSges(7,2) - mrSges(2,1);
t133 = m(6) * pkin(11) - mrSges(5,2) + mrSges(6,3);
t56 = sin(qJ(4));
t61 = cos(qJ(4));
t132 = -pkin(9) * m(6) - mrSges(6,1) * t61 + mrSges(6,2) * t56 - mrSges(5,1);
t131 = t56 * mrSges(6,1) + t61 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(7,3) + mrSges(8,3) + mrSges(9,3);
t130 = 0.1e1 / pkin(2);
t127 = -pkin(8) - pkin(10);
t126 = -pkin(8) + pkin(10);
t125 = pkin(1) * t63;
t118 = -pkin(8) ^ 2 + pkin(10) ^ 2;
t117 = cos(pkin(19));
t116 = sin(pkin(19));
t115 = cos(pkin(17));
t114 = sin(pkin(17));
t113 = pkin(18) + pkin(19);
t53 = pkin(12) + r_base(3);
t59 = sin(qJ(1));
t112 = t59 * pkin(13) + r_base(2);
t64 = cos(qJ(1));
t111 = t64 * pkin(13) + r_base(1);
t108 = pkin(1) - t124;
t41 = t59 * t125 + t112;
t42 = t64 * t125 + t111;
t47 = t58 * pkin(1) + t53;
t107 = cos(t113);
t106 = sin(t113);
t102 = t110 - t119;
t96 = t130 * (-pkin(5) * t121 + t102 * t108);
t94 = -t96 / 0.2e1;
t95 = t130 * (pkin(5) * t102 * t46 + t108 * t29) / 0.2e1;
t21 = (t116 * t95 + t117 * t94) * t31;
t93 = t96 / 0.2e1;
t22 = (t116 * t93 + t117 * t95) * t31;
t20 = t21 * t63 - t22 * t58;
t19 = t21 * t58 + t22 * t63;
t57 = sin(qJ(3));
t62 = cos(qJ(3));
t100 = t57 * t63 + t58 * t62;
t43 = t57 * t58 - t62 * t63;
t92 = t57 * t93 + t62 * t95;
t91 = t57 * t95 + t62 * t94;
t90 = t31 * (-t106 * t92 - t107 * t91);
t89 = t31 * (t106 * t91 - t107 * t92);
t88 = pkin(4) * t90;
t87 = pkin(3) - t88;
t86 = -pkin(3) * t90 + pkin(4);
t85 = -0.2e1 * pkin(3) * t88 + pkin(4) ^ 2;
t84 = pkin(3) ^ 2 + t85;
t83 = 0.1e1 / t84;
t82 = 0.1e1 / pkin(8) * t83;
t81 = 0.1e1 / pkin(10) * t83;
t80 = t84 + t118;
t79 = t84 - t118;
t78 = sqrt(-((pkin(3) - t126) * (pkin(3) + t126) + t85) * ((pkin(3) - t127) * (pkin(3) + t127) + t85));
t77 = t78 * t89;
t76 = (-pkin(4) * t77 + t79 * t87) * t82;
t75 = (-pkin(3) * t77 + t80 * t86) * t81;
t74 = -(pkin(4) * t79 * t89 + t78 * t87) * t82 / 0.2e1;
t73 = (pkin(3) * t80 * t89 + t78 * t86) * t81 / 0.2e1;
t55 = cos(pkin(18));
t54 = sin(pkin(18));
t37 = t43 * t64;
t36 = t100 * t64;
t35 = t43 * t59;
t34 = t100 * t59;
t18 = t19 * t64;
t17 = t20 * t64;
t16 = t19 * t59;
t15 = t20 * t59;
t10 = -t55 * t76 / 0.2e1 + t54 * t74;
t9 = t54 * t76 / 0.2e1 + t55 * t74;
t8 = t115 * t73 + t114 * t75 / 0.2e1;
t7 = -t115 * t75 / 0.2e1 + t114 * t73;
t1 = (-m(7) * pkin(14) - mrSges(3,1) * t58 - mrSges(3,2) * t63 + mrSges(4,1) * t100 - (t19 * t55 - t20 * t54) * t139 - m(1) * r_base(3) - mrSges(8,1) * t19 - mrSges(8,2) * t20 - t23 * mrSges(7,2) - t24 * mrSges(7,1) - (t10 * t19 + t20 * t9) * mrSges(9,1) - (t10 * t20 - t19 * t9) * mrSges(9,2) - mrSges(4,2) * t43 - mrSges(1,3) - mrSges(2,3) + t137 * (-pkin(4) * t100 + t47) + t133 * (t100 * t8 + t43 * t7) + (-m(3) + t138) * t53 + t132 * (-t100 * t7 + t43 * t8) + t135 * t47) * g(3) + (-m(3) * t112 - mrSges(4,1) * t35 - mrSges(8,1) * t15 + mrSges(8,2) * t16 - mrSges(4,2) * t34 - (t15 * t55 + t16 * t54) * t139 - (t10 * t15 - t16 * t9) * mrSges(9,1) - (-t10 * t16 - t15 * t9) * mrSges(9,2) - mrSges(1,2) + t135 * t41 + t137 * (t35 * pkin(4) + t41) + t133 * (t34 * t7 - t35 * t8) + t136 * r_base(2) + t134 * t59 + t132 * (t34 * t8 + t35 * t7) + t131 * t64) * g(2) + (-m(3) * t111 - mrSges(4,1) * t37 - (t17 * t55 + t18 * t54) * t139 - mrSges(8,1) * t17 + mrSges(8,2) * t18 - (t10 * t17 - t18 * t9) * mrSges(9,1) - (-t10 * t18 - t17 * t9) * mrSges(9,2) - mrSges(4,2) * t36 - mrSges(1,1) + t135 * t42 + t137 * (t37 * pkin(4) + t42) + t133 * (t36 * t7 - t37 * t8) + t136 * r_base(1) + t134 * t64 + t132 * (t36 * t8 + t37 * t7) - t131 * t59) * g(1);
U = t1;
