% Calculate joint inertia matrix for
% palh1m1OL
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
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1OL_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_inertiaJ_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_inertiaJ_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_inertiaJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_inertiaJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_inertiaJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:40
% EndTime: 2020-04-15 19:28:47
% DurationCPUTime: 1.12s
% Computational Cost: add. (1066->233), mult. (2127->354), div. (0->0), fcn. (2311->20), ass. (0->115)
t107 = cos(qJ(5));
t99 = sin(qJ(5));
t122 = mrSges(6,1) * t107 - mrSges(6,2) * t99;
t152 = pkin(9) * t122;
t90 = t99 ^ 2;
t148 = mrSges(6,3) * t90;
t84 = pkin(11) * t148;
t91 = t107 ^ 2;
t147 = mrSges(6,3) * t91;
t85 = pkin(11) * t147;
t158 = t152 + t84 + t85;
t108 = cos(qJ(4));
t150 = pkin(5) * t108;
t83 = -pkin(9) - t150;
t48 = t83 * t122;
t100 = sin(qJ(4));
t81 = pkin(5) * t100 + pkin(11);
t67 = t81 * t148;
t68 = t81 * t147;
t86 = mrSges(5,1) * t150;
t157 = -t48 + t67 + t68 + t86;
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t109 = cos(qJ(3));
t110 = cos(qJ(2));
t65 = -t101 * t102 + t109 * t110;
t78 = t102 * pkin(1) - pkin(15);
t39 = -pkin(5) * t65 + t78;
t155 = 0.2e1 * t39;
t151 = pkin(1) * t109;
t95 = sin(qJ(9));
t149 = mrSges(10,2) * t95;
t146 = Ifges(6,4) * t99;
t62 = t101 * t110 + t102 * t109;
t34 = t100 * t65 + t108 * t62;
t145 = t34 * t99;
t82 = pkin(1) * t101 + pkin(5);
t45 = t100 * t82 - t108 * t151;
t144 = t45 * mrSges(5,2);
t138 = sin(qJ(10));
t92 = sin(pkin(19));
t93 = cos(pkin(19));
t94 = cos(qJ(10));
t56 = -t138 * t93 + t92 * t94;
t57 = -t138 * t92 - t93 * t94;
t105 = cos(qJ(7));
t97 = sin(qJ(7));
t61 = -t102 * t105 - t110 * t97;
t64 = -t102 * t97 + t105 * t110;
t12 = -t56 * t64 + t57 * t61;
t13 = t56 * t61 + t57 * t64;
t143 = Ifges(11,5) * t13 + Ifges(11,6) * t12;
t134 = t107 * t34;
t32 = t100 * t62 - t108 * t65;
t142 = Ifges(6,5) * t134 + Ifges(6,3) * t32;
t103 = cos(qJ(9));
t104 = cos(qJ(8));
t96 = sin(qJ(8));
t60 = -t102 * t104 - t110 * t96;
t63 = -t102 * t96 + t104 * t110;
t31 = -t103 * t60 + t63 * t95;
t33 = -t103 * t63 - t60 * t95;
t141 = Ifges(10,5) * t33 + Ifges(10,6) * t31;
t140 = Ifges(6,5) * t99 + Ifges(6,6) * t107;
t139 = t90 + t91;
t137 = mrSges(10,1) * t103;
t136 = Ifges(6,4) * t107;
t135 = t100 * mrSges(5,2);
t69 = pkin(1) * t97 + pkin(4) * t92;
t71 = pkin(1) * t105 + pkin(4) * t93;
t16 = t56 * t71 + t57 * t69;
t133 = t16 * mrSges(11,2);
t24 = (t56 * t93 + t57 * t92) * pkin(4);
t132 = t24 * mrSges(11,2);
t131 = Ifges(11,3) + Ifges(8,3);
t23 = (-t56 * t92 + t57 * t93) * pkin(4);
t21 = t23 * mrSges(11,1);
t130 = Ifges(11,3) + t21;
t129 = pkin(5) * t135;
t72 = Ifges(6,2) * t107 + t146;
t73 = Ifges(6,1) * t99 + t136;
t128 = t107 * t72 + t99 * t73 + Ifges(5,3);
t127 = t139 * t81;
t126 = Ifges(4,3) + t128;
t125 = Ifges(8,5) * t64 + Ifges(8,6) * t61 + t143;
t36 = Ifges(10,3) + Ifges(9,3) + (0.2e1 * t149 - 0.2e1 * t137 + m(10) * (t103 ^ 2 + t95 ^ 2) * pkin(2)) * pkin(2);
t8 = -mrSges(6,2) * t32 - mrSges(6,3) * t145;
t9 = mrSges(6,1) * t32 - mrSges(6,3) * t134;
t124 = t107 * t8 - t9 * t99;
t123 = mrSges(6,1) * t99 + mrSges(6,2) * t107;
t121 = mrSges(8,1) * t105 - mrSges(8,2) * t97;
t120 = mrSges(4,1) * t101 + mrSges(4,2) * t109;
t3 = Ifges(6,6) * t32 + (-Ifges(6,2) * t99 + t136) * t34;
t4 = Ifges(6,5) * t32 + (Ifges(6,1) * t107 - t146) * t34;
t119 = t99 * t4 / 0.2e1 - t72 * t145 / 0.2e1 + t73 * t134 / 0.2e1 + t107 * t3 / 0.2e1 + Ifges(5,5) * t34 + (t140 / 0.2e1 - Ifges(5,6)) * t32;
t44 = t100 * t151 + t108 * t82;
t118 = Ifges(4,5) * t62 + Ifges(4,6) * t65 + t119;
t117 = Ifges(9,5) * t63 + Ifges(9,6) * t60 + t141 + (t103 * t33 - t31 * t95) * pkin(2) * mrSges(10,3);
t42 = -pkin(9) - t44;
t35 = t42 * t122;
t43 = pkin(11) + t45;
t37 = t43 * t148;
t38 = t43 * t147;
t40 = t44 * mrSges(5,1);
t116 = t128 - t35 + t37 + t38 + t40 - t144;
t106 = cos(qJ(6));
t98 = sin(qJ(6));
t66 = Ifges(10,3) + (-t137 + t149) * pkin(2);
t46 = -pkin(2) * t60 - pkin(15);
t20 = (-t61 * t93 - t64 * t92) * pkin(4) + t78;
t15 = -t56 * t69 + t57 * t71;
t14 = t15 * mrSges(11,1);
t7 = t123 * t34;
t6 = pkin(9) * t32 - pkin(11) * t34 + t39;
t1 = [0.2e1 * (-mrSges(4,1) * t65 - mrSges(8,1) * t61 + mrSges(4,2) * t62 + mrSges(8,2) * t64) * t78 - 0.2e1 * (mrSges(3,1) * t102 - mrSges(9,1) * t60 + mrSges(3,2) * t110 + mrSges(9,2) * t63) * pkin(15) + (m(6) * t139 * t6 + 0.2e1 * t107 * t9 + 0.2e1 * t8 * t99) * t6 + (m(8) + m(4)) * t78 ^ 2 + (m(9) + m(3)) * pkin(15) ^ 2 + Ifges(9,1) * t63 ^ 2 + Ifges(8,1) * t64 ^ 2 + Ifges(4,2) * t65 ^ 2 + Ifges(3,1) * t110 ^ 2 + (0.2e1 * Ifges(9,4) * t63 + Ifges(9,2) * t60) * t60 + t106 * (Ifges(7,4) * t98 + Ifges(7,2) * t106) + t98 * (Ifges(7,1) * t98 + Ifges(7,4) * t106) + 0.2e1 * pkin(14) * (-mrSges(7,1) * t106 + mrSges(7,2) * t98) + (-0.2e1 * Ifges(3,4) * t110 + Ifges(3,2) * t102) * t102 + m(5) * t39 ^ 2 + 0.2e1 * t46 * (-mrSges(10,1) * t31 + mrSges(10,2) * t33) + m(10) * t46 ^ 2 + t33 * (Ifges(10,1) * t33 + Ifges(10,4) * t31) + t31 * (Ifges(10,4) * t33 + Ifges(10,2) * t31) + 0.2e1 * t20 * (-mrSges(11,1) * t12 + mrSges(11,2) * t13) + m(11) * t20 ^ 2 + t13 * (Ifges(11,1) * t13 + Ifges(11,4) * t12) + t12 * (Ifges(11,4) * t13 + Ifges(11,2) * t12) + (Ifges(4,1) * t62 + 0.2e1 * Ifges(4,4) * t65) * t62 + (0.2e1 * Ifges(8,4) * t64 + Ifges(8,2) * t61) * t61 + (mrSges(5,2) * t155 + Ifges(5,1) * t34 + t107 * t4 - t3 * t99 + (-Ifges(6,6) * t99 - (2 * Ifges(5,4))) * t32) * t34 + (mrSges(5,1) * t155 + Ifges(5,2) * t32 + t142) * t32 + Ifges(2,3) + m(7) * pkin(14) ^ 2; ((-t105 * t64 + t61 * t97) * mrSges(8,3) + (-t101 * t62 - t109 * t65) * mrSges(4,3)) * pkin(1) + t124 * t43 + (-t32 * t45 - t34 * t44) * mrSges(5,3) + (t12 * t16 - t13 * t15) * mrSges(11,3) + Ifges(3,5) * t110 - Ifges(3,6) * t102 + t42 * t7 + t118 + t117 + t125; m(6) * (t139 * t43 ^ 2 + t42 ^ 2) + t126 + m(5) * (t44 ^ 2 + t45 ^ 2) + m(11) * (t15 ^ 2 + t16 ^ 2) + 0.2e1 * t40 + 0.2e1 * t38 + t36 + 0.2e1 * t37 - 0.2e1 * t35 + 0.2e1 * t14 - 0.2e1 * t144 - 0.2e1 * t133 + Ifges(3,3) + t131 + (0.2e1 * t120 + 0.2e1 * t121 + (m(4) * (t101 ^ 2 + t109 ^ 2) + m(8) * (t105 ^ 2 + t97 ^ 2)) * pkin(1)) * pkin(1); t7 * t83 + t124 * t81 + (-t100 * t32 - t108 * t34) * mrSges(5,3) * pkin(5) + t118; m(6) * (t127 * t43 + t42 * t83) + (m(5) * (t100 * t45 + t108 * t44) - t135) * pkin(5) + t120 * pkin(1) + t116 + Ifges(4,3) + t157; -0.2e1 * t129 - 0.2e1 * t48 + 0.2e1 * t67 + 0.2e1 * t68 + 0.2e1 * t86 + m(6) * (t139 * t81 ^ 2 + t83 ^ 2) + m(5) * (t100 ^ 2 + t108 ^ 2) * pkin(5) ^ 2 + t126; -pkin(9) * t7 + pkin(11) * t124 + t119; m(6) * (pkin(11) * t139 * t43 - pkin(9) * t42) + t116 + t158; -t129 + m(6) * (-pkin(9) * t83 + pkin(11) * t127) + t128 + t157 + t158; 0.2e1 * t84 + 0.2e1 * t85 + 0.2e1 * t152 + m(6) * (pkin(11) ^ 2 * t139 + pkin(9) ^ 2) + t128; -Ifges(6,6) * t145 + t122 * t6 + t142; -t123 * t43 + t140; -t123 * t81 + t140; -pkin(11) * t123 + t140; Ifges(6,3); Ifges(7,5) * t98 + Ifges(7,6) * t106; 0; 0; 0; 0; Ifges(7,3); (t12 * t24 - t13 * t23) * mrSges(11,3) + t125; Ifges(8,3) + m(11) * (t15 * t23 + t16 * t24) + t14 + (-t24 - t16) * mrSges(11,2) + t121 * pkin(1) + t130; 0; 0; 0; 0; m(11) * (t23 ^ 2 + t24 ^ 2) - 0.2e1 * t132 + 0.2e1 * t21 + t131; t117; t36; 0; 0; 0; 0; 0; t36; t141; t66; 0; 0; 0; 0; 0; t66; Ifges(10,3); t143; Ifges(11,3) + t14 - t133; 0; 0; 0; 0; t130 - t132; 0; 0; Ifges(11,3); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_13_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46), t1(56), t1(67), t1(79); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47), t1(57), t1(68), t1(80); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48), t1(58), t1(69), t1(81); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49), t1(59), t1(70), t1(82); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50), t1(60), t1(71), t1(83); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51), t1(61), t1(72), t1(84); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52), t1(62), t1(73), t1(85); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53), t1(63), t1(74), t1(86); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54), t1(64), t1(75), t1(87); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55), t1(65), t1(76), t1(88); t1(56), t1(57), t1(58), t1(59), t1(60), t1(61), t1(62), t1(63), t1(64), t1(65), t1(66), t1(77), t1(89); t1(67), t1(68), t1(69), t1(70), t1(71), t1(72), t1(73), t1(74), t1(75), t1(76), t1(77), t1(78), t1(90); t1(79), t1(80), t1(81), t1(82), t1(83), t1(84), t1(85), t1(86), t1(87), t1(88), t1(89), t1(90), t1(91);];
Mq = res;
