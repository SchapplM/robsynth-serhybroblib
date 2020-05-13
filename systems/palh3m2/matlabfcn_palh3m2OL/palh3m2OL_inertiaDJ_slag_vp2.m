% Calculate time derivative of joint inertia matrix for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
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
% MqD [10x10]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh3m2OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),zeros(16,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_inertiaDJ_slag_vp2: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_inertiaDJ_slag_vp2: qJD has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_inertiaDJ_slag_vp2: pkin has to be [16x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2OL_inertiaDJ_slag_vp2: m has to be [9x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [9,3]), ...
  'palh3m2OL_inertiaDJ_slag_vp2: mrSges has to be [9x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [9 6]), ...
  'palh3m2OL_inertiaDJ_slag_vp2: Ifges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:33:44
% EndTime: 2020-05-07 04:34:01
% DurationCPUTime: 2.56s
% Computational Cost: add. (2911->288), mult. (6932->474), div. (0->0), fcn. (6855->16), ass. (0->152)
t113 = sin(qJ(5));
t107 = t113 ^ 2;
t120 = cos(qJ(5));
t108 = t120 ^ 2;
t156 = t107 + t108;
t111 = sin(qJ(7));
t116 = sin(qJ(2));
t118 = cos(qJ(7));
t123 = cos(qJ(2));
t83 = -t111 * t116 + t118 * t123;
t196 = -0.2e1 * t83;
t153 = qJD(5) * t120;
t114 = sin(qJ(4));
t121 = cos(qJ(4));
t115 = sin(qJ(3));
t122 = cos(qJ(3));
t134 = t115 * t123 + t116 * t122;
t84 = t115 * t116 - t122 * t123;
t57 = -t114 * t134 - t121 * t84;
t193 = qJD(2) + qJD(3);
t66 = t193 * t84;
t67 = t193 * t134;
t24 = -qJD(4) * t57 + t114 * t67 + t121 * t66;
t58 = t114 * t84 - t121 * t134;
t133 = t113 * t24 + t58 * t153;
t106 = qJD(2) * t116 * pkin(1);
t50 = -pkin(4) * t67 + t106;
t195 = 0.2e1 * t50;
t102 = -pkin(1) * t123 - pkin(12);
t186 = 0.2e1 * t102;
t109 = sin(pkin(15));
t110 = cos(pkin(15));
t117 = cos(qJ(8));
t182 = sin(qJ(8));
t82 = -t109 * t182 - t110 * t117;
t194 = t156 * t121;
t192 = qJD(2) + qJD(7);
t104 = -pkin(4) * t121 - pkin(8);
t139 = mrSges(6,1) * t113 + mrSges(6,2) * t120;
t87 = t139 * qJD(5);
t69 = t104 * t87;
t173 = pkin(4) * qJD(4);
t172 = mrSges(6,1) * t120;
t91 = mrSges(6,2) * t113 - t172;
t79 = t114 * t91 * t173;
t155 = qJD(4) * t121;
t150 = pkin(4) * t155;
t141 = mrSges(6,3) * t150;
t96 = t107 * t141;
t97 = t108 * t141;
t191 = t69 + t79 + t96 + t97;
t25 = qJD(4) * t58 + t114 * t66 - t121 * t67;
t8 = pkin(8) * t25 - pkin(10) * t24 + t50;
t190 = 0.2e1 * t8;
t85 = t111 * t123 + t116 * t118;
t188 = 0.2e1 * t85;
t187 = -0.2e1 * t134;
t184 = pkin(8) * t87;
t130 = -t109 * t117 + t110 * t182;
t64 = t192 * t83;
t65 = t192 * t85;
t75 = t82 * qJD(8);
t76 = t130 * qJD(8);
t11 = t130 * t65 + t64 * t82 + t75 * t83 + t76 * t85;
t12 = t130 * t64 - t65 * t82 - t75 * t85 + t76 * t83;
t183 = Ifges(9,5) * t11 + Ifges(9,6) * t12;
t181 = Ifges(6,5) * t25;
t180 = Ifges(6,5) * t57;
t179 = Ifges(6,6) * t57;
t178 = t25 * Ifges(6,6);
t174 = pkin(1) * qJD(7);
t94 = pkin(1) * t111 + pkin(3) * t109;
t95 = pkin(1) * t118 + pkin(3) * t110;
t27 = t75 * t95 + t76 * t94 + (t111 * t130 + t118 * t82) * t174;
t177 = t27 * mrSges(9,2);
t158 = t114 * t115;
t52 = t150 + (qJD(3) + qJD(4)) * pkin(1) * (-t121 * t122 + t158);
t176 = t52 * mrSges(5,2);
t164 = t120 * t24;
t175 = Ifges(6,5) * t164 + Ifges(6,3) * t25;
t171 = Ifges(6,4) * t113;
t170 = Ifges(6,4) * t120;
t169 = Ifges(6,6) * t113;
t168 = t107 * t52;
t167 = t108 * t52;
t165 = t113 * t58;
t163 = t120 * t58;
t157 = t115 * t121;
t154 = qJD(5) * t113;
t68 = -pkin(4) * t84 + t102;
t29 = pkin(8) * t57 - pkin(10) * t58 + t68;
t152 = 0.2e1 * t29;
t151 = 0.2e1 * t123;
t148 = pkin(1) * t122 - pkin(4);
t46 = (-t109 * t75 + t110 * t76) * pkin(3);
t44 = t46 * mrSges(9,1);
t45 = (t109 * t76 + t110 * t75) * pkin(3);
t146 = -t45 * mrSges(9,2) + t44;
t137 = -Ifges(6,2) * t113 + t170;
t18 = t137 * t58 + t179;
t145 = -t18 - t179;
t143 = Ifges(8,5) * t64 - Ifges(8,6) * t65 + t183;
t142 = t156 * t52;
t98 = t114 * t148;
t73 = -pkin(1) * t157 - t98;
t140 = -t114 * mrSges(5,1) - t121 * mrSges(5,2);
t138 = Ifges(6,1) * t120 - t171;
t32 = -mrSges(6,2) * t57 - mrSges(6,3) * t165;
t33 = mrSges(6,1) * t57 - mrSges(6,3) * t163;
t136 = -t113 * t33 + t120 * t32;
t132 = t154 * t58 - t164;
t88 = t137 * qJD(5);
t89 = t138 * qJD(5);
t92 = Ifges(6,2) * t120 + t171;
t93 = Ifges(6,1) * t113 + t170;
t131 = t113 * t89 + t120 * t88 + t153 * t93 - t154 * t92;
t129 = qJD(3) * (mrSges(4,1) * t115 + mrSges(4,2) * t122);
t128 = qJD(7) * (-mrSges(8,1) * t111 - mrSges(8,2) * t118);
t72 = pkin(1) * t158 - t121 * t148;
t6 = mrSges(6,1) * t25 + mrSges(6,3) * t132;
t7 = -mrSges(6,2) * t25 - mrSges(6,3) * t133;
t127 = -t113 * t6 + t120 * t7 + (-t113 * t32 - t120 * t33) * qJD(5);
t105 = Ifges(6,5) * t153;
t19 = t138 * t58 + t180;
t3 = -Ifges(6,4) * t132 - Ifges(6,2) * t133 + t178;
t4 = -Ifges(6,1) * t132 - Ifges(6,4) * t133 + t181;
t126 = t113 * t4 / 0.2e1 + t93 * t164 / 0.2e1 + t19 * t153 / 0.2e1 + t25 * (Ifges(6,5) * t113 + Ifges(6,6) * t120) / 0.2e1 + t120 * t3 / 0.2e1 + Ifges(5,5) * t24 - t88 * t165 / 0.2e1 + t89 * t163 / 0.2e1 + t57 * (-Ifges(6,6) * t154 + t105) / 0.2e1 - Ifges(5,6) * t25 - t133 * t92 / 0.2e1 - (t58 * t93 + t18) * t154 / 0.2e1;
t53 = qJD(4) * t98 + (t115 * t155 + (t114 * t122 + t157) * qJD(3)) * pkin(1);
t43 = t53 * t91;
t47 = mrSges(6,3) * t168;
t48 = mrSges(6,3) * t167;
t51 = t53 * mrSges(5,1);
t70 = -pkin(8) - t72;
t63 = t70 * t87;
t125 = t131 - t43 + t47 + t48 + t51 + t63 - t176;
t124 = Ifges(4,5) * t66 + Ifges(4,6) * t67 + t126;
t119 = cos(qJ(6));
t112 = sin(qJ(6));
t103 = pkin(4) * t114 + pkin(10);
t71 = pkin(10) + t73;
t55 = (t109 * t82 - t110 * t130) * pkin(3);
t54 = (t109 * t130 + t110 * t82) * pkin(3);
t42 = -t130 * t95 + t82 * t94;
t41 = t130 * t94 + t82 * t95;
t35 = -t130 * t83 + t82 * t85;
t34 = t130 * t85 + t82 * t83;
t31 = t139 * t58;
t30 = t106 + (-t109 * t64 + t110 * t65) * pkin(3);
t28 = -t75 * t94 + t76 * t95 + (-t111 * t82 + t118 * t130) * t174;
t26 = t28 * mrSges(9,1);
t5 = mrSges(6,1) * t133 - mrSges(6,2) * t132;
t1 = [0.2e1 * t84 * Ifges(4,2) * t67 + 0.2e1 * t58 * t24 * Ifges(5,1) + t66 * Ifges(4,1) * t187 + 0.2e1 * t34 * Ifges(9,2) * t12 + 0.2e1 * t35 * t11 * Ifges(9,1) + 0.2e1 * t30 * (-mrSges(9,1) * t34 + mrSges(9,2) * t35) + 0.2e1 * t68 * (mrSges(5,1) * t25 + mrSges(5,2) * t24) + Ifges(8,2) * t65 * t196 + t64 * Ifges(8,1) * t188 + (mrSges(8,1) * t65 + mrSges(8,2) * t64) * t186 + (-mrSges(4,1) * t67 + mrSges(4,2) * t66) * t186 + (m(5) * t68 + mrSges(5,2) * t58) * t195 + (mrSges(5,1) * t195 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t25 + t175) * t57 + m(6) * t156 * t29 * t190 + 0.2e1 * (t11 * t34 + t12 * t35) * Ifges(9,4) + 0.2e1 * (t64 * t83 - t65 * t85) * Ifges(8,4) + 0.2e1 * (-t24 * t57 - t25 * t58) * Ifges(5,4) + 0.2e1 * (-t134 * t67 + t66 * t84) * Ifges(4,4) + (t19 * t24 + t33 * t190 + (qJD(5) * t32 + t6) * t152 + (qJD(5) * t145 + t181 + t4) * t58) * t120 + (t32 * t190 + (-qJD(5) * t33 + t7) * t152 + t145 * t24 + (-t178 - t3 + (-t19 - t180) * qJD(5)) * t58) * t113 + 0.2e1 * (pkin(6) * (mrSges(7,1) * t112 + mrSges(7,2) * t119) + (-t112 ^ 2 + t119 ^ 2) * Ifges(7,4) + (-Ifges(7,2) + Ifges(7,1)) * t112 * t119) * qJD(6) + ((-pkin(12) * mrSges(3,2) + Ifges(3,4) * t123) * t151 + (-0.2e1 * pkin(12) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t116 + (Ifges(3,1) - Ifges(3,2)) * t151 + (-0.2e1 * mrSges(4,1) * t84 + mrSges(8,1) * t196 + mrSges(4,2) * t187 + mrSges(8,2) * t188 + (m(4) + m(8)) * t186) * pkin(1)) * t116) * qJD(2) + 0.2e1 * (m(9) * t30 - mrSges(9,1) * t12 + mrSges(9,2) * t11) * ((-t109 * t85 - t110 * t83) * pkin(3) + t102); t124 + ((-t111 * t65 - t118 * t64 + (t111 * t85 + t118 * t83) * qJD(7)) * mrSges(8,3) + (-t115 * t67 + t122 * t66 + (t115 * t134 - t122 * t84) * qJD(3)) * mrSges(4,3)) * pkin(1) + t127 * t71 + t136 * t52 + (Ifges(3,5) * t123 - Ifges(3,6) * t116) * qJD(2) + (-t11 * t41 + t12 * t42 + t27 * t34 - t28 * t35) * mrSges(9,3) + (-t24 * t72 - t25 * t73 - t52 * t57 - t53 * t58) * mrSges(5,3) - t53 * t31 + t70 * t5 + t143; -0.2e1 * t176 - 0.2e1 * t177 + 0.2e1 * t26 - 0.2e1 * t43 + 0.2e1 * t47 + 0.2e1 * t48 + 0.2e1 * t51 + 0.2e1 * t63 + 0.2e1 * m(9) * (t27 * t42 + t28 * t41) + 0.2e1 * m(6) * (t142 * t71 - t53 * t70) + 0.2e1 * m(5) * (t52 * t73 + t53 * t72) + t131 + 0.2e1 * (t128 + t129) * pkin(1); t124 + ((-t114 * t25 - t121 * t24) * mrSges(5,3) + ((mrSges(5,3) * t58 + t31) * t114 + (-mrSges(5,3) * t57 + t136) * t121) * qJD(4)) * pkin(4) + t127 * t103 + t104 * t5; m(6) * (-t104 * t53 + (t167 + t168) * t103) + (m(5) * (t114 * t52 + t121 * t53) + (m(6) * (t114 * t70 + t194 * t71) + m(5) * (-t114 * t72 + t121 * t73) + t140) * qJD(4)) * pkin(4) + t125 + pkin(1) * t129 + t191; 0.2e1 * t69 + 0.2e1 * t79 + 0.2e1 * t96 + 0.2e1 * t97 + 0.2e1 * (m(6) * (t103 * t194 + t104 * t114) + t140) * t173 + t131; -pkin(8) * t5 + pkin(10) * t127 + t126; -t184 + m(6) * (pkin(8) * t53 + pkin(10) * t142) + t125; -t184 + (m(6) * (-pkin(8) * t114 + pkin(10) * t194) + t140) * t173 + t131 + t191; t131 - 0.2e1 * t184; t8 * t172 + (-mrSges(6,2) * t8 - Ifges(6,6) * t24) * t113 + ((-mrSges(6,2) * t29 - Ifges(6,6) * t58) * t120 + (-mrSges(6,1) * t29 - Ifges(6,5) * t58) * t113) * qJD(5) + t175; t105 - t139 * t52 + (-t71 * t172 + (mrSges(6,2) * t71 - Ifges(6,6)) * t113) * qJD(5); t105 - t139 * t150 + (t103 * t91 - t169) * qJD(5); t105 + (pkin(10) * t91 - t169) * qJD(5); 0; (Ifges(7,5) * t119 - Ifges(7,6) * t112) * qJD(6); 0; 0; 0; 0; 0; (-t11 * t54 + t12 * t55 + t34 * t45 - t35 * t46) * mrSges(9,3) + t143; t26 + m(9) * (t27 * t55 + t28 * t54 + t41 * t46 + t42 * t45) + t44 + (-t27 - t45) * mrSges(9,2) + pkin(1) * t128; 0; 0; 0; 0; 0.2e1 * m(9) * (t45 * t55 + t46 * t54) + 0.2e1 * t146; t183; t26 - t177; 0; 0; 0; 0; t146; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_10_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55);];
Mq = res;
