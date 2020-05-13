% Calculate Gravitation load on the joints for
% picker2Dm2IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 09:21
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = picker2Dm2IC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm2IC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [11x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 09:21:01
% EndTime: 2020-05-11 09:21:02
% DurationCPUTime: 0.98s
% Computational Cost: add. (998->187), mult. (467->182), div. (37->10), fcn. (360->59), ass. (0->106)
t78 = qJ(1) + qJ(2);
t73 = qJ(3) + t78;
t54 = sin(t73);
t57 = cos(t73);
t109 = t57 * mrSges(4,1) - t54 * mrSges(4,2);
t114 = qJ(9) + t78;
t61 = qJ(3) + t114;
t40 = sin(t61);
t43 = cos(t61);
t110 = -t43 * mrSges(10,1) + t40 * mrSges(10,2);
t33 = pkin(6) * t57;
t131 = -m(10) * t33 - t109 - t110;
t84 = cos(qJ(1));
t130 = pkin(1) * t84;
t117 = qJ(2) + qJ(4);
t72 = qJ(1) + t117;
t53 = sin(t72);
t34 = pkin(4) * t53;
t85 = 0.1e1 / pkin(5);
t116 = qJ(3) - qJ(6);
t75 = (pkin(8) + qJ(5));
t103 = t75 + t116;
t77 = qJ(1) + qJ(8);
t91 = t103 - t77;
t104 = t75 - t116;
t92 = t104 - t77;
t129 = 0.1e1 / (pkin(6) * (cos(t92) - cos(t91)) + (-cos(qJ(9) - t92) + cos(qJ(9) + t91)) * pkin(2)) * t85;
t113 = -qJ(9) - t116;
t8 = 0.1e1 / ((-cos(qJ(4) - t116) + cos(qJ(4) + t116)) * pkin(6) + (cos(qJ(4) + t113) - cos(qJ(4) - t113)) * pkin(2));
t87 = 0.1e1 / pkin(3);
t128 = t8 * t87;
t82 = sin(qJ(1));
t74 = t82 * pkin(1);
t127 = t57 * mrSges(4,2);
t60 = -qJ(9) + t72;
t25 = t60 - t116;
t59 = qJ(4) + t114;
t26 = t59 + t116;
t126 = sin(t26) - sin(t25);
t125 = cos(t26) - cos(t25);
t124 = t40 * mrSges(10,1) + t43 * mrSges(10,2);
t56 = cos(t72);
t123 = t53 * mrSges(5,1) + t56 * mrSges(5,2);
t122 = -cos(0.2e1 * t77) + cos((2 * t75));
t121 = sin(t60) - sin(t59);
t120 = cos(t60) - cos(t59);
t66 = sin(t78);
t68 = cos(t78);
t119 = -t66 * mrSges(3,1) - t68 * mrSges(3,2);
t50 = pkin(3) * t66;
t118 = t50 + t74;
t80 = sin(qJ(4));
t76 = 0.1e1 / t80;
t86 = 0.1e1 / pkin(4);
t115 = pkin(1) * t76 * t86;
t51 = pkin(2) * t66;
t112 = -pkin(6) * t54 + t51;
t111 = -qJ(1) + t75;
t71 = qJ(6) + t78;
t52 = sin(t71);
t55 = cos(t71);
t14 = t55 * mrSges(7,1) - t52 * mrSges(7,2);
t65 = sin(t77);
t67 = cos(t77);
t108 = t67 * mrSges(9,1) - t65 * mrSges(9,2);
t58 = qJ(10) + t72;
t35 = sin(t58);
t36 = cos(t58);
t107 = -t36 * mrSges(11,1) + t35 * mrSges(11,2);
t106 = -qJ(4) + t111;
t105 = qJ(4) + t111;
t102 = -pkin(2) * t68 - t130;
t101 = -pkin(3) * t68 - t130;
t100 = -t54 * mrSges(4,1) - t127;
t13 = -t52 * mrSges(7,1) - t55 * mrSges(7,2);
t99 = -t65 * mrSges(9,1) - t67 * mrSges(9,2);
t98 = t53 * t68 - t56 * t66;
t79 = sin(qJ(7));
t83 = cos(qJ(7));
t97 = t53 * t79 + t56 * t83;
t95 = t35 * mrSges(11,1) + t36 * mrSges(11,2);
t94 = -qJ(8) + t106;
t93 = -qJ(8) + t105;
t89 = t95 - t123;
t27 = t53 * mrSges(5,2);
t88 = -t27 + (pkin(4) * m(11) + mrSges(5,1)) * t56 + t107;
t81 = sin(qJ(2));
t70 = -qJ(9) + t75;
t69 = qJ(9) + t75;
t64 = sin(t117);
t63 = cos(t75);
t62 = sin(t75);
t46 = t66 * mrSges(3,2);
t45 = qJ(9) + t103;
t44 = -qJ(9) + t104;
t16 = -pkin(5) * t65 + t74;
t15 = -pkin(5) * t67 + t130;
t12 = pkin(4) * t56 - t101;
t11 = t34 + t118;
t7 = -g(1) * t13 - g(2) * t14;
t6 = -g(1) * t124 - g(2) * t110;
t5 = g(1) * t95 + g(2) * t107;
t3 = t131 * g(2) + (t127 + (m(10) * pkin(6) + mrSges(4,1)) * t54 - t124) * g(1);
t2 = t88 * g(2) + (-m(11) * t34 + t89) * g(1);
t1 = (-t14 - t46 + (mrSges(3,1) + (m(11) + m(5)) * pkin(3) + (m(4) + m(10)) * pkin(2)) * t68 + t88 + t131) * g(2) + (-m(10) * t112 - m(11) * (t50 + t34) - t13 + (-pkin(2) * m(4) - pkin(3) * m(5)) * t66 + t89 - t100 + t119 - t124) * g(1);
t4 = [(-t122 * pkin(5) + (-cos(qJ(8) + 0.2e1 * qJ(1)) + cos((2 * pkin(8)) + (2 * qJ(5)) + qJ(8))) * pkin(1)) * t85 / t122 * (-g(1) * t99 - g(2) * t108) - g(2) * (m(5) * t101 - t56 * mrSges(5,1) + t27) - g(1) * (m(7) * t74 + t13) - g(2) * (-m(7) * t130 + t14) - g(1) * (m(9) * t74 + t99) - g(1) * (t82 * mrSges(2,1) + t84 * mrSges(2,2)) - g(2) * (-t84 * mrSges(2,1) + t82 * mrSges(2,2)) - g(1) * (m(3) * t74 - t119) - g(2) * (-m(9) * t130 + t108) - g(1) * (m(10) * (t112 + t74) + t124) - g(2) * (m(10) * (t102 + t33) + t110) - g(1) * (m(11) * t11 - t95) - g(2) * (-m(11) * t12 - t107) - g(2) * (-m(3) * t130 - t68 * mrSges(3,1) + t46) - g(1) * (m(4) * (t51 + t74) + t100) - g(2) * (m(4) * t102 + t109) - g(1) * (m(5) * t118 + t123) - t81 * t5 * t115 + sin(qJ(8)) / sin(-qJ(8) + t111) * (-g(1) * (-mrSges(6,1) * t62 - mrSges(6,2) * t63) - g(2) * (mrSges(6,1) * t63 - mrSges(6,2) * t62)) + (((t121 * t11 + t120 * t12) * t128 + (-(sin(t70) - sin(t69)) * t16 - (cos(t70) - cos(t69)) * t15) * t129) * t7 + ((-t126 * t11 - t125 * t12) * t128 + ((sin(t45) - sin(t44)) * t16 + (cos(t45) - cos(t44)) * t15) * t129) * t3) * pkin(2) + ((-pkin(1) * t64 - pkin(3) * t80) * t76 * t1 - (pkin(3) * (cos(t106) - cos(t105)) + (cos(-qJ(2) + t94) - cos(qJ(2) + t93)) * pkin(5)) * pkin(1) * t85 / (cos(t94) - cos(t93)) * t6 + (pkin(3) * t81 + pkin(4) * t64) * t2 * t115) * t87; -g(1) * (mrSges(8,1) * t83 - mrSges(8,2) * t79) - g(2) * (mrSges(8,1) * t79 + mrSges(8,2) * t83) + ((-t121 * t7 + t126 * t3) * t83 + (t120 * t7 - t125 * t3) * t79) * t8 * pkin(2) + (((t97 * t2 + t98 * t5) * pkin(4) + (t2 - t5) * pkin(3) * (t66 * t79 + t68 * t83)) * t86 - (t1 + t6) * t97) / t98;];
taug = t4(:);
