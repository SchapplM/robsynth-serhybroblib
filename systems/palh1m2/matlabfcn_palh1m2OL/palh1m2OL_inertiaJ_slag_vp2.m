% Calculate joint inertia matrix for
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
% m [11x1]
%   mass of all robot links (including the base)
% mrSges [11x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [13x13]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_inertiaJ_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_inertiaJ_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:37
% EndTime: 2020-05-02 21:17:38
% DurationCPUTime: 1.08s
% Computational Cost: add. (1066->236), mult. (2127->358), div. (0->0), fcn. (2303->20), ass. (0->117)
t100 = sin(qJ(5));
t108 = cos(qJ(5));
t141 = mrSges(6,1) * t108;
t69 = mrSges(6,2) * t100 - t141;
t155 = pkin(9) * t69;
t91 = t100 ^ 2;
t151 = mrSges(6,3) * t91;
t85 = pkin(11) * t151;
t92 = t108 ^ 2;
t150 = mrSges(6,3) * t92;
t86 = pkin(11) * t150;
t160 = t85 + t86 - t155;
t109 = cos(qJ(4));
t153 = pkin(5) * t109;
t84 = -pkin(9) - t153;
t48 = t84 * t69;
t101 = sin(qJ(4));
t83 = pkin(5) * t101 + pkin(11);
t67 = t83 * t151;
t68 = t83 * t150;
t87 = mrSges(5,1) * t153;
t159 = t48 + t67 + t68 + t87;
t110 = cos(qJ(3));
t154 = pkin(1) * t110;
t96 = sin(qJ(9));
t152 = mrSges(10,2) * t96;
t102 = sin(qJ(3));
t103 = sin(qJ(2));
t111 = cos(qJ(2));
t64 = t102 * t111 + t103 * t110;
t130 = t110 * t111;
t65 = -t102 * t103 + t130;
t33 = -t101 * t64 + t109 * t65;
t149 = Ifges(6,6) * t33;
t82 = pkin(1) * t102 + pkin(5);
t44 = t82 * t101 - t109 * t154;
t148 = t44 * mrSges(5,2);
t142 = sin(qJ(10));
t93 = sin(pkin(19));
t94 = cos(pkin(19));
t95 = cos(qJ(10));
t56 = -t142 * t94 + t93 * t95;
t57 = -t142 * t93 - t94 * t95;
t106 = cos(qJ(7));
t98 = sin(qJ(7));
t61 = -t103 * t98 + t106 * t111;
t62 = -t103 * t106 - t111 * t98;
t12 = -t56 * t61 + t57 * t62;
t13 = t56 * t62 + t57 * t61;
t147 = Ifges(11,5) * t13 + Ifges(11,6) * t12;
t34 = t101 * t65 + t109 * t64;
t135 = t108 * t34;
t146 = Ifges(6,5) * t135 - Ifges(6,3) * t33;
t104 = cos(qJ(9));
t105 = cos(qJ(8));
t97 = sin(qJ(8));
t60 = -t103 * t105 - t111 * t97;
t63 = -t103 * t97 + t105 * t111;
t31 = -t104 * t60 + t63 * t96;
t32 = -t104 * t63 - t60 * t96;
t145 = Ifges(10,5) * t32 + Ifges(10,6) * t31;
t144 = Ifges(6,5) * t100 + Ifges(6,6) * t108;
t79 = t103 * pkin(1) - pkin(15);
t143 = t91 + t92;
t140 = mrSges(10,1) * t104;
t139 = Ifges(6,4) * t100;
t138 = Ifges(6,4) * t108;
t137 = t100 * t34;
t136 = t101 * mrSges(5,2);
t72 = pkin(1) * t98 + pkin(4) * t93;
t73 = pkin(1) * t106 + pkin(4) * t94;
t16 = t56 * t73 + t57 * t72;
t134 = t16 * mrSges(11,2);
t24 = (t56 * t94 + t57 * t93) * pkin(4);
t133 = t24 * mrSges(11,2);
t132 = Ifges(11,3) + Ifges(8,3);
t23 = (-t56 * t93 + t57 * t94) * pkin(4);
t21 = t23 * mrSges(11,1);
t131 = Ifges(11,3) + t21;
t129 = pkin(5) * t136;
t70 = Ifges(6,2) * t108 + t139;
t71 = Ifges(6,1) * t100 + t138;
t128 = t100 * t71 + t108 * t70 + Ifges(5,3);
t127 = t143 * t83;
t126 = Ifges(4,3) + t128;
t125 = Ifges(8,5) * t61 + Ifges(8,6) * t62 + t147;
t36 = Ifges(10,3) + Ifges(9,3) + (0.2e1 * t152 - 0.2e1 * t140 + m(10) * (t104 ^ 2 + t96 ^ 2) * pkin(2)) * pkin(2);
t8 = mrSges(6,2) * t33 - mrSges(6,3) * t137;
t9 = -mrSges(6,1) * t33 - mrSges(6,3) * t135;
t124 = -t100 * t9 + t108 * t8;
t123 = mrSges(8,1) * t106 - mrSges(8,2) * t98;
t122 = mrSges(4,1) * t102 + mrSges(4,2) * t110;
t121 = mrSges(6,1) * t100 + mrSges(6,2) * t108;
t41 = -pkin(5) * t130 - pkin(15) + (pkin(5) * t102 + pkin(1)) * t103;
t3 = -t149 + (-Ifges(6,2) * t100 + t138) * t34;
t4 = -Ifges(6,5) * t33 + (Ifges(6,1) * t108 - t139) * t34;
t120 = t100 * t4 / 0.2e1 - t70 * t137 / 0.2e1 + t71 * t135 / 0.2e1 + t108 * t3 / 0.2e1 + Ifges(5,5) * t34 + (-t144 / 0.2e1 + Ifges(5,6)) * t33;
t45 = t101 * t154 + t109 * t82;
t119 = Ifges(4,5) * t64 + Ifges(4,6) * t65 + t120;
t118 = Ifges(9,5) * t63 + Ifges(9,6) * t60 + t145 + (t104 * t32 - t31 * t96) * pkin(2) * mrSges(10,3);
t43 = -pkin(9) - t45;
t35 = t43 * t69;
t42 = pkin(11) + t44;
t37 = t42 * t151;
t38 = t42 * t150;
t39 = t45 * mrSges(5,1);
t117 = t128 + t35 + t37 + t38 + t39 - t148;
t107 = cos(qJ(6));
t99 = sin(qJ(6));
t66 = Ifges(10,3) + (-t140 + t152) * pkin(2);
t46 = -pkin(2) * t60 - pkin(15);
t20 = (-t61 * t93 - t62 * t94) * pkin(4) + t79;
t15 = -t56 * t72 + t57 * t73;
t14 = t15 * mrSges(11,1);
t7 = t121 * t34;
t6 = -pkin(9) * t33 - pkin(11) * t34 + t41;
t1 = [(m(8) + m(4)) * t79 ^ 2 + (m(9) + m(3)) * pkin(15) ^ 2 + Ifges(9,1) * t63 ^ 2 + Ifges(8,2) * t62 ^ 2 + Ifges(4,2) * t65 ^ 2 + Ifges(3,1) * t111 ^ 2 + (0.2e1 * Ifges(9,4) * t63 + Ifges(9,2) * t60) * t60 + t99 * (Ifges(7,1) * t99 + Ifges(7,4) * t107) + t107 * (Ifges(7,4) * t99 + Ifges(7,2) * t107) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t107 + mrSges(7,2) * t99) + 0.2e1 * (-mrSges(4,1) * t65 - mrSges(8,1) * t62 + mrSges(4,2) * t64 + mrSges(8,2) * t61) * t79 - 0.2e1 * (mrSges(3,1) * t103 - mrSges(9,1) * t60 + mrSges(3,2) * t111 + mrSges(9,2) * t63) * pkin(15) + (m(6) * t143 * t6 + 0.2e1 * t100 * t8 + 0.2e1 * t108 * t9) * t6 + (-0.2e1 * Ifges(3,4) * t111 + Ifges(3,2) * t103) * t103 + 0.2e1 * t46 * (-mrSges(10,1) * t31 + mrSges(10,2) * t32) + m(10) * t46 ^ 2 + m(5) * t41 ^ 2 + 0.2e1 * t20 * (-mrSges(11,1) * t12 + mrSges(11,2) * t13) + m(11) * t20 ^ 2 + t32 * (Ifges(10,1) * t32 + Ifges(10,4) * t31) + t31 * (Ifges(10,4) * t32 + Ifges(10,2) * t31) + t12 * (Ifges(11,4) * t13 + Ifges(11,2) * t12) + t13 * (Ifges(11,1) * t13 + Ifges(11,4) * t12) + (Ifges(4,1) * t64 + 0.2e1 * Ifges(4,4) * t65) * t64 + (Ifges(8,1) * t61 + 0.2e1 * Ifges(8,4) * t62) * t61 + (-0.2e1 * mrSges(5,1) * t41 + Ifges(5,2) * t33 - t146) * t33 + (0.2e1 * mrSges(5,2) * t41 + Ifges(5,1) * t34 + 0.2e1 * Ifges(5,4) * t33 + t108 * t4 + (-t3 + t149) * t100) * t34 + m(7) * pkin(14) ^ 2 + Ifges(2,3); ((-t106 * t61 + t62 * t98) * mrSges(8,3) + (-t102 * t64 - t110 * t65) * mrSges(4,3)) * pkin(1) + t124 * t42 + (t33 * t44 - t34 * t45) * mrSges(5,3) + (t12 * t16 - t13 * t15) * mrSges(11,3) + Ifges(3,5) * t111 - Ifges(3,6) * t103 + t43 * t7 + t119 + t118 + t125; t126 + m(6) * (t143 * t42 ^ 2 + t43 ^ 2) + m(5) * (t44 ^ 2 + t45 ^ 2) + m(11) * (t15 ^ 2 + t16 ^ 2) + 0.2e1 * t38 + 0.2e1 * t39 + 0.2e1 * t37 + 0.2e1 * t35 + t36 + 0.2e1 * t14 - 0.2e1 * t148 - 0.2e1 * t134 + Ifges(3,3) + t132 + (0.2e1 * t122 + 0.2e1 * t123 + (m(4) * (t102 ^ 2 + t110 ^ 2) + m(8) * (t106 ^ 2 + t98 ^ 2)) * pkin(1)) * pkin(1); t7 * t84 + t124 * t83 + (t101 * t33 - t109 * t34) * mrSges(5,3) * pkin(5) + t119; m(6) * (t127 * t42 + t43 * t84) + (m(5) * (t101 * t44 + t109 * t45) - t136) * pkin(5) + t122 * pkin(1) + t117 + Ifges(4,3) + t159; -0.2e1 * t129 + 0.2e1 * t48 + 0.2e1 * t67 + 0.2e1 * t68 + 0.2e1 * t87 + m(6) * (t143 * t83 ^ 2 + t84 ^ 2) + m(5) * (t101 ^ 2 + t109 ^ 2) * pkin(5) ^ 2 + t126; -pkin(9) * t7 + pkin(11) * t124 + t120; m(6) * (pkin(11) * t143 * t42 - pkin(9) * t43) + t117 + t160; -t129 + m(6) * (-pkin(9) * t84 + pkin(11) * t127) + t128 + t159 + t160; 0.2e1 * t86 + 0.2e1 * t85 - 0.2e1 * t155 + m(6) * (pkin(11) ^ 2 * t143 + pkin(9) ^ 2) + t128; t6 * t141 + (-mrSges(6,2) * t6 - Ifges(6,6) * t34) * t100 + t146; -t121 * t42 + t144; -t121 * t83 + t144; -pkin(11) * t121 + t144; Ifges(6,3); Ifges(7,5) * t99 + Ifges(7,6) * t107; 0; 0; 0; 0; Ifges(7,3); (t12 * t24 - t13 * t23) * mrSges(11,3) + t125; Ifges(8,3) + t14 + m(11) * (t15 * t23 + t16 * t24) + (-t16 - t24) * mrSges(11,2) + t123 * pkin(1) + t131; 0; 0; 0; 0; m(11) * (t23 ^ 2 + t24 ^ 2) - 0.2e1 * t133 + 0.2e1 * t21 + t132; t118; t36; 0; 0; 0; 0; 0; t36; t145; t66; 0; 0; 0; 0; 0; t66; Ifges(10,3); t147; Ifges(11,3) + t14 - t134; 0; 0; 0; 0; t131 - t133; 0; 0; Ifges(11,3); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_13_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46), t1(56), t1(67), t1(79); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47), t1(57), t1(68), t1(80); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48), t1(58), t1(69), t1(81); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49), t1(59), t1(70), t1(82); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50), t1(60), t1(71), t1(83); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51), t1(61), t1(72), t1(84); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52), t1(62), t1(73), t1(85); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53), t1(63), t1(74), t1(86); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54), t1(64), t1(75), t1(87); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55), t1(65), t1(76), t1(88); t1(56), t1(57), t1(58), t1(59), t1(60), t1(61), t1(62), t1(63), t1(64), t1(65), t1(66), t1(77), t1(89); t1(67), t1(68), t1(69), t1(70), t1(71), t1(72), t1(73), t1(74), t1(75), t1(76), t1(77), t1(78), t1(90); t1(79), t1(80), t1(81), t1(82), t1(83), t1(84), t1(85), t1(86), t1(87), t1(88), t1(89), t1(90), t1(91);];
Mq = res;
