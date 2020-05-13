% Calculate joint inertia matrix for
% palh3m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% mrSges [9x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [9x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [10x10]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:16
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m1OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1OL_inertiaJ_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1OL_inertiaJ_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m1OL_inertiaJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m1OL_inertiaJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m1OL_inertiaJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:04:02
% EndTime: 2020-04-20 17:04:09
% DurationCPUTime: 0.83s
% Computational Cost: add. (950->202), mult. (1887->307), div. (0->0), fcn. (2055->16), ass. (0->100)
t80 = sin(qJ(5));
t87 = cos(qJ(5));
t99 = mrSges(6,1) * t87 - mrSges(6,2) * t80;
t124 = pkin(8) * t99;
t74 = t80 ^ 2;
t122 = mrSges(6,3) * t74;
t69 = pkin(10) * t122;
t75 = t87 ^ 2;
t121 = mrSges(6,3) * t75;
t70 = pkin(10) * t121;
t130 = t124 + t69 + t70;
t88 = cos(qJ(4));
t125 = pkin(4) * t88;
t67 = -pkin(8) - t125;
t41 = t67 * t99;
t81 = sin(qJ(4));
t66 = pkin(4) * t81 + pkin(10);
t54 = t66 * t122;
t55 = t66 * t121;
t71 = mrSges(5,1) * t125;
t129 = -t41 + t54 + t55 + t71;
t82 = sin(qJ(3));
t83 = sin(qJ(2));
t89 = cos(qJ(3));
t90 = cos(qJ(2));
t51 = t82 * t83 - t89 * t90;
t65 = -pkin(1) * t90 - pkin(12);
t33 = -pkin(4) * t51 + t65;
t128 = 0.2e1 * t33;
t126 = pkin(1) * t82;
t123 = sin(qJ(8));
t120 = Ifges(6,4) * t80;
t119 = Ifges(6,4) * t87;
t76 = sin(pkin(15));
t77 = cos(pkin(15));
t84 = cos(qJ(8));
t46 = -t123 * t77 + t76 * t84;
t47 = -t123 * t76 - t77 * t84;
t78 = sin(qJ(7));
t59 = pkin(1) * t78 + pkin(3) * t76;
t85 = cos(qJ(7));
t60 = pkin(1) * t85 + pkin(3) * t77;
t18 = t46 * t60 + t47 * t59;
t118 = t18 * mrSges(9,2);
t24 = (t46 * t77 + t47 * t76) * pkin(3);
t117 = t24 * mrSges(9,2);
t53 = -t82 * t90 - t83 * t89;
t29 = t51 * t81 + t53 * t88;
t116 = t29 * t80;
t115 = t29 * t87;
t68 = -pkin(1) * t89 + pkin(4);
t39 = -t126 * t88 + t81 * t68;
t114 = t39 * mrSges(5,2);
t113 = t81 * mrSges(5,2);
t112 = Ifges(8,3) + Ifges(9,3);
t23 = (-t46 * t76 + t47 * t77) * pkin(3);
t21 = t23 * mrSges(9,1);
t111 = Ifges(9,3) + t21;
t50 = -t78 * t83 + t85 * t90;
t52 = t78 * t90 + t83 * t85;
t12 = -t46 * t52 + t47 * t50;
t13 = t46 * t50 + t47 * t52;
t110 = Ifges(9,5) * t13 + Ifges(9,6) * t12;
t28 = -t88 * t51 + t53 * t81;
t109 = Ifges(6,5) * t115 + Ifges(6,3) * t28;
t108 = Ifges(6,5) * t80 + Ifges(6,6) * t87;
t107 = t74 + t75;
t106 = pkin(4) * t113;
t57 = Ifges(6,2) * t87 + t120;
t58 = Ifges(6,1) * t80 + t119;
t105 = t87 * t57 + t80 * t58 + Ifges(5,3);
t104 = t107 * t66;
t103 = Ifges(4,3) + t105;
t102 = Ifges(8,5) * t52 + Ifges(8,6) * t50 + t110;
t8 = -mrSges(6,2) * t28 - mrSges(6,3) * t116;
t9 = mrSges(6,1) * t28 - mrSges(6,3) * t115;
t101 = t8 * t87 - t80 * t9;
t100 = -mrSges(4,1) * t89 + mrSges(4,2) * t82;
t98 = mrSges(6,1) * t80 + mrSges(6,2) * t87;
t97 = mrSges(8,1) * t85 - mrSges(8,2) * t78;
t38 = t126 * t81 + t68 * t88;
t3 = Ifges(6,6) * t28 + (-Ifges(6,2) * t80 + t119) * t29;
t4 = Ifges(6,5) * t28 + (Ifges(6,1) * t87 - t120) * t29;
t96 = t80 * t4 / 0.2e1 - t57 * t116 / 0.2e1 + t58 * t115 / 0.2e1 + t87 * t3 / 0.2e1 + Ifges(5,5) * t29 + (t108 / 0.2e1 - Ifges(5,6)) * t28;
t95 = Ifges(4,5) * t53 + Ifges(4,6) * t51 + t96;
t36 = -pkin(8) - t38;
t30 = t36 * t99;
t37 = pkin(10) + t39;
t31 = t37 * t122;
t32 = t37 * t121;
t34 = t38 * mrSges(5,1);
t94 = t105 - t30 + t31 + t32 + t34 - t114;
t86 = cos(qJ(6));
t79 = sin(qJ(6));
t20 = (-t50 * t77 - t52 * t76) * pkin(3) + t65;
t17 = -t46 * t59 + t47 * t60;
t14 = t17 * mrSges(9,1);
t7 = t98 * t29;
t6 = pkin(8) * t28 - pkin(10) * t29 + t33;
t1 = [-0.2e1 * pkin(12) * (-mrSges(3,1) * t90 + mrSges(3,2) * t83) + t83 * (Ifges(3,1) * t83 + Ifges(3,4) * t90) + t90 * (Ifges(3,4) * t83 + Ifges(3,2) * t90) + t86 * (Ifges(7,4) * t79 + Ifges(7,2) * t86) + t79 * (Ifges(7,1) * t79 + Ifges(7,4) * t86) + 0.2e1 * pkin(6) * (-mrSges(7,1) * t86 + mrSges(7,2) * t79) + Ifges(4,1) * t53 ^ 2 + Ifges(8,1) * t52 ^ 2 + m(5) * t33 ^ 2 + 0.2e1 * t20 * (-mrSges(9,1) * t12 + mrSges(9,2) * t13) + m(9) * t20 ^ 2 + t13 * (Ifges(9,1) * t13 + Ifges(9,4) * t12) + t12 * (Ifges(9,4) * t13 + Ifges(9,2) * t12) + (mrSges(5,2) * t128 + Ifges(5,1) * t29 - t3 * t80 + t4 * t87 + (-Ifges(6,6) * t80 - (2 * Ifges(5,4))) * t28) * t29 + (mrSges(5,1) * t128 + Ifges(5,2) * t28 + t109) * t28 + (m(8) + m(4)) * t65 ^ 2 + m(7) * pkin(6) ^ 2 + m(3) * pkin(12) ^ 2 + Ifges(2,3) + (0.2e1 * Ifges(4,4) * t53 + Ifges(4,2) * t51) * t51 + (0.2e1 * Ifges(8,4) * t52 + Ifges(8,2) * t50) * t50 + (m(6) * t107 * t6 + 0.2e1 * t8 * t80 + 0.2e1 * t87 * t9) * t6 + 0.2e1 * (-mrSges(4,1) * t51 - mrSges(8,1) * t50 + mrSges(4,2) * t53 + mrSges(8,2) * t52) * t65; ((t50 * t78 - t52 * t85) * mrSges(8,3) + (-t51 * t82 + t53 * t89) * mrSges(4,3)) * pkin(1) + (t12 * t18 - t13 * t17) * mrSges(9,3) + t101 * t37 + (-t28 * t39 - t29 * t38) * mrSges(5,3) + Ifges(3,6) * t90 + Ifges(3,5) * t83 + t36 * t7 + t95 + t102; -0.2e1 * t114 - 0.2e1 * t118 + Ifges(3,3) + 0.2e1 * t14 - 0.2e1 * t30 + 0.2e1 * t31 + 0.2e1 * t32 + 0.2e1 * t34 + m(6) * (t107 * t37 ^ 2 + t36 ^ 2) + m(5) * (t38 ^ 2 + t39 ^ 2) + m(9) * (t17 ^ 2 + t18 ^ 2) + t103 + t112 + (0.2e1 * t100 + 0.2e1 * t97 + (m(8) * (t78 ^ 2 + t85 ^ 2) + m(4) * (t82 ^ 2 + t89 ^ 2)) * pkin(1)) * pkin(1); t67 * t7 + t101 * t66 + (-t28 * t81 - t29 * t88) * mrSges(5,3) * pkin(4) + t95; t100 * pkin(1) + t94 + m(6) * (t104 * t37 + t36 * t67) + (m(5) * (t38 * t88 + t39 * t81) - t113) * pkin(4) + Ifges(4,3) + t129; -0.2e1 * t106 - 0.2e1 * t41 + 0.2e1 * t54 + 0.2e1 * t55 + 0.2e1 * t71 + m(6) * (t107 * t66 ^ 2 + t67 ^ 2) + m(5) * (t81 ^ 2 + t88 ^ 2) * pkin(4) ^ 2 + t103; -pkin(8) * t7 + pkin(10) * t101 + t96; m(6) * (pkin(10) * t107 * t37 - pkin(8) * t36) + t94 + t130; m(6) * (-pkin(8) * t67 + pkin(10) * t104) - t106 + t105 + t129 + t130; m(6) * (pkin(10) ^ 2 * t107 + pkin(8) ^ 2) + 0.2e1 * t124 + 0.2e1 * t69 + 0.2e1 * t70 + t105; -Ifges(6,6) * t116 + t6 * t99 + t109; -t37 * t98 + t108; -t66 * t98 + t108; -pkin(10) * t98 + t108; Ifges(6,3); Ifges(7,5) * t79 + Ifges(7,6) * t86; 0; 0; 0; 0; Ifges(7,3); (t12 * t24 - t13 * t23) * mrSges(9,3) + t102; t14 + m(9) * (t17 * t23 + t18 * t24) + Ifges(8,3) + (-t24 - t18) * mrSges(9,2) + t97 * pkin(1) + t111; 0; 0; 0; 0; -0.2e1 * t117 + 0.2e1 * t21 + m(9) * (t23 ^ 2 + t24 ^ 2) + t112; t110; Ifges(9,3) + t14 - t118; 0; 0; 0; 0; t111 - t117; Ifges(9,3); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_10_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55);];
Mq = res;
