% Calculate time derivative of joint inertia matrix for
% palh1m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% qJD [13x1]
%   Generalized joint velocities
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
% MqD [13x13]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-15 19:46
% Revision: 2d0abd6fcc3afe6f578a07ad3d897ec57baa6ba1 (2020-04-13)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m1OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m1OL_inertiaDJ_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m1OL_inertiaDJ_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m1OL_inertiaDJ_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m1OL_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m1OL_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m1OL_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-15 19:28:41
% EndTime: 2020-04-15 19:29:00
% DurationCPUTime: 3.52s
% Computational Cost: add. (3242->322), mult. (7717->535), div. (0->0), fcn. (7641->20), ass. (0->173)
t134 = sin(qJ(5));
t125 = t134 ^ 2;
t142 = cos(qJ(5));
t126 = t142 ^ 2;
t182 = t125 + t126;
t136 = sin(qJ(3));
t137 = sin(qJ(2));
t144 = cos(qJ(3));
t145 = cos(qJ(2));
t102 = -t136 * t137 + t144 * t145;
t228 = -0.2e1 * t102;
t180 = qJD(5) * t142;
t135 = sin(qJ(4));
t143 = cos(qJ(4));
t99 = t136 * t145 + t144 * t137;
t160 = t143 * t102 - t135 * t99;
t224 = qJD(2) + qJD(3);
t77 = t224 * t99;
t78 = t224 * t102;
t29 = qJD(4) * t160 - t135 * t77 + t143 * t78;
t65 = t102 * t135 + t143 * t99;
t156 = t134 * t29 + t65 * t180;
t123 = qJD(2) * t145 * pkin(1);
t55 = pkin(5) * t77 + t123;
t227 = 0.2e1 * t55;
t130 = sin(qJ(9));
t138 = cos(qJ(9));
t94 = (mrSges(10,1) * t130 + mrSges(10,2) * t138) * qJD(9) * pkin(2);
t87 = 0.2e1 * t94;
t119 = t137 * pkin(1) - pkin(15);
t214 = 0.2e1 * t119;
t127 = sin(pkin(19));
t128 = cos(pkin(19));
t129 = cos(qJ(10));
t201 = sin(qJ(10));
t96 = -t127 * t201 - t128 * t129;
t225 = t182 * t143;
t223 = qJD(2) + qJD(7);
t222 = qJD(2) + qJD(8);
t221 = qJD(3) + qJD(4);
t199 = pkin(5) * qJD(4);
t176 = t143 * t199;
t167 = mrSges(6,3) * t176;
t112 = t125 * t167;
t113 = t126 * t167;
t165 = mrSges(6,1) * t134 + mrSges(6,2) * t142;
t103 = t165 * qJD(5);
t121 = -pkin(5) * t143 - pkin(9);
t80 = t121 * t103;
t198 = mrSges(6,1) * t142;
t108 = mrSges(6,2) * t134 - t198;
t92 = t135 * t108 * t199;
t220 = t112 + t113 + t80 + t92;
t30 = qJD(4) * t65 + t135 * t78 + t143 * t77;
t8 = pkin(9) * t30 - pkin(11) * t29 + t55;
t219 = 0.2e1 * t8;
t218 = -0.2e1 * pkin(15);
t216 = 0.2e1 * t99;
t132 = sin(qJ(7));
t140 = cos(qJ(7));
t158 = t132 * t137 - t140 * t145;
t215 = -0.2e1 * t158;
t153 = -t127 * t129 + t128 * t201;
t98 = -t132 * t145 - t140 * t137;
t75 = t223 * t98;
t76 = t223 * t158;
t83 = t96 * qJD(10);
t84 = t153 * qJD(10);
t11 = -t153 * t76 - t158 * t84 + t75 * t96 + t83 * t98;
t12 = t153 * t75 + t158 * t83 + t76 * t96 + t84 * t98;
t211 = Ifges(11,5) * t11 + Ifges(11,6) * t12;
t210 = pkin(1) * t136;
t209 = pkin(9) * t103;
t208 = Ifges(6,5) * t30;
t207 = Ifges(6,5) * t160;
t206 = Ifges(6,6) * t160;
t205 = t30 * Ifges(6,6);
t183 = t135 * t144;
t59 = t176 + t221 * pkin(1) * (t136 * t143 + t183);
t204 = t59 * mrSges(5,2);
t190 = t142 * t29;
t203 = Ifges(6,5) * t190 + Ifges(6,3) * t30;
t131 = sin(qJ(8));
t139 = cos(qJ(8));
t159 = t131 * t137 - t139 * t145;
t97 = -t131 * t145 - t139 * t137;
t62 = -t130 * t159 - t138 * t97;
t73 = t222 * t97;
t74 = t222 * t159;
t27 = qJD(9) * t62 - t130 * t74 - t138 * t73;
t161 = t130 * t97 - t138 * t159;
t28 = qJD(9) * t161 + t130 * t73 - t138 * t74;
t202 = Ifges(10,5) * t27 + Ifges(10,6) * t28;
t200 = pkin(1) * qJD(7);
t197 = Ifges(6,4) * t134;
t196 = Ifges(6,4) * t142;
t195 = Ifges(6,6) * t134;
t194 = t125 * t59;
t193 = t126 * t59;
t191 = t134 * t65;
t189 = t142 * t65;
t107 = pkin(1) * t132 + pkin(4) * t127;
t109 = pkin(1) * t140 + pkin(4) * t128;
t32 = t107 * t84 + t109 * t83 + (t132 * t153 + t140 * t96) * t200;
t188 = t32 * mrSges(11,2);
t181 = qJD(5) * t134;
t79 = -pkin(5) * t102 + t119;
t34 = -pkin(9) * t160 - pkin(11) * t65 + t79;
t179 = 0.2e1 * t34;
t178 = 0.2e1 * t137;
t177 = pkin(1) * t143 * t144;
t174 = pkin(5) + t210;
t163 = -Ifges(6,2) * t134 + t196;
t18 = t163 * t65 - t206;
t173 = -t18 + t206;
t171 = -t181 / 0.2e1;
t170 = Ifges(8,5) * t75 + Ifges(8,6) * t76 + t211;
t51 = (-t127 * t83 + t128 * t84) * pkin(4);
t49 = t51 * mrSges(11,1);
t50 = (t127 * t84 + t128 * t83) * pkin(4);
t169 = -t50 * mrSges(11,2) + t49;
t168 = t182 * t59;
t114 = t135 * t174;
t86 = t114 - t177;
t166 = -t135 * mrSges(5,1) - t143 * mrSges(5,2);
t164 = Ifges(6,1) * t142 - t197;
t37 = mrSges(6,2) * t160 - mrSges(6,3) * t191;
t38 = -mrSges(6,1) * t160 - mrSges(6,3) * t189;
t162 = -t134 * t38 + t142 * t37;
t155 = t181 * t65 - t190;
t104 = t163 * qJD(5);
t105 = t164 * qJD(5);
t110 = Ifges(6,2) * t142 + t197;
t111 = Ifges(6,1) * t134 + t196;
t154 = t142 * t104 + t134 * t105 - t110 * t181 + t111 * t180;
t152 = qJD(3) * (mrSges(4,1) * t144 - mrSges(4,2) * t136);
t151 = qJD(7) * (-mrSges(8,1) * t132 - mrSges(8,2) * t140);
t85 = pkin(1) * t183 + t143 * t174;
t6 = mrSges(6,1) * t30 + mrSges(6,3) * t155;
t7 = -mrSges(6,2) * t30 - mrSges(6,3) * t156;
t150 = -t134 * t6 + t142 * t7 + (-t134 * t37 - t142 * t38) * qJD(5);
t122 = Ifges(6,5) * t180;
t19 = t164 * t65 - t207;
t3 = -Ifges(6,4) * t155 - Ifges(6,2) * t156 + t205;
t4 = -Ifges(6,1) * t155 - Ifges(6,4) * t156 + t208;
t149 = t134 * t4 / 0.2e1 + t18 * t171 + t19 * t180 / 0.2e1 + t30 * (Ifges(6,5) * t134 + Ifges(6,6) * t142) / 0.2e1 + t142 * t3 / 0.2e1 + Ifges(5,5) * t29 - t104 * t191 / 0.2e1 + t105 * t189 / 0.2e1 - t160 * (-Ifges(6,6) * t181 + t122) / 0.2e1 - Ifges(5,6) * t30 - t156 * t110 / 0.2e1 + (t190 / 0.2e1 + t65 * t171) * t111;
t60 = -t135 * qJD(3) * t210 - qJD(4) * t114 + t177 * t221;
t48 = t60 * t108;
t53 = mrSges(6,3) * t194;
t54 = mrSges(6,3) * t193;
t56 = t60 * mrSges(5,1);
t81 = -pkin(9) - t85;
t72 = t81 * t103;
t148 = t154 - t48 + t53 + t54 + t56 + t72 - t204;
t147 = Ifges(4,5) * t78 - Ifges(4,6) * t77 + t149;
t146 = Ifges(9,5) * t73 + Ifges(9,6) * t74 + t202 + (t138 * t27 - t130 * t28 + (t130 * t161 - t138 * t62) * qJD(9)) * pkin(2) * mrSges(10,3);
t141 = cos(qJ(6));
t133 = sin(qJ(6));
t120 = pkin(5) * t135 + pkin(11);
t88 = -pkin(2) * t97 - pkin(15);
t82 = pkin(11) + t86;
t58 = (t127 * t96 - t128 * t153) * pkin(4);
t57 = (t127 * t153 + t128 * t96) * pkin(4);
t44 = t107 * t96 - t109 * t153;
t43 = t107 * t153 + t109 * t96;
t40 = -t153 * t98 - t158 * t96;
t39 = -t153 * t158 + t96 * t98;
t36 = t165 * t65;
t35 = t123 + (-t127 * t75 - t128 * t76) * pkin(4);
t33 = -t107 * t83 + t109 * t84 + (-t132 * t96 + t140 * t153) * t200;
t31 = t33 * mrSges(11,1);
t5 = mrSges(6,1) * t156 - mrSges(6,2) * t155;
t1 = [((pkin(15) * mrSges(3,2) + Ifges(3,4) * t137) * t178 + (mrSges(3,1) * t218 - 0.2e1 * Ifges(3,4) * t145 + (Ifges(3,2) - Ifges(3,1)) * t178 + (mrSges(4,1) * t228 - 0.2e1 * mrSges(8,1) * t98 + mrSges(4,2) * t216 + mrSges(8,2) * t215 + (m(4) + m(8)) * t214) * pkin(1)) * t145) * qJD(2) + (-mrSges(8,1) * t76 + mrSges(8,2) * t75) * t214 + (-mrSges(9,1) * t74 + mrSges(9,2) * t73) * t218 + 0.2e1 * (pkin(14) * (mrSges(7,1) * t133 + mrSges(7,2) * t141) + (-t133 ^ 2 + t141 ^ 2) * Ifges(7,4) + (-Ifges(7,2) + Ifges(7,1)) * t133 * t141) * qJD(6) + 0.2e1 * (t11 * t39 + t12 * t40) * Ifges(11,4) + 0.2e1 * (m(11) * t35 - mrSges(11,1) * t12 + mrSges(11,2) * t11) * ((t127 * t158 - t128 * t98) * pkin(4) + t119) + 0.2e1 * (t160 * t29 - t30 * t65) * Ifges(5,4) - (mrSges(5,1) * t227 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t30 + t203) * t160 + 0.2e1 * ((-m(10) * t88 + mrSges(10,1) * t62 + mrSges(10,2) * t161) * pkin(2) + Ifges(9,2) * t97) * t74 + 0.2e1 * (-t161 * t28 + t62 * t27) * Ifges(10,4) + 0.2e1 * (-t159 * t74 + t73 * t97) * Ifges(9,4) + 0.2e1 * (-t158 * t76 + t75 * t98) * Ifges(8,4) + (mrSges(4,1) * t77 + mrSges(4,2) * t78) * t214 + 0.2e1 * (t102 * t78 - t77 * t99) * Ifges(4,4) + 0.2e1 * t62 * Ifges(10,2) * t28 + (m(5) * t79 + mrSges(5,2) * t65) * t227 + 0.2e1 * t98 * Ifges(8,2) * t76 + 0.2e1 * t39 * Ifges(11,2) * t12 + 0.2e1 * t29 * t65 * Ifges(5,1) + 0.2e1 * t11 * t40 * Ifges(11,1) + (t37 * t219 + (-qJD(5) * t38 + t7) * t179 + t173 * t29 + (-t205 - t3 + (-t19 + t207) * qJD(5)) * t65) * t134 + 0.2e1 * t79 * (mrSges(5,1) * t30 + mrSges(5,2) * t29) + 0.2e1 * t88 * (-mrSges(10,1) * t28 + mrSges(10,2) * t27) + 0.2e1 * t35 * (-mrSges(11,1) * t39 + mrSges(11,2) * t40) + m(6) * t182 * t34 * t219 - 0.2e1 * t27 * t161 * Ifges(10,1) - 0.2e1 * t73 * t159 * Ifges(9,1) + t75 * Ifges(8,1) * t215 + t78 * Ifges(4,1) * t216 + t77 * Ifges(4,2) * t228 + (t19 * t29 + t38 * t219 + (qJD(5) * t37 + t6) * t179 + (qJD(5) * t173 + t208 + t4) * t65) * t142; t147 + ((t132 * t76 - t140 * t75 + (-t132 * t158 + t140 * t98) * qJD(7)) * mrSges(8,3) + (-t136 * t78 + t144 * t77 + (t102 * t136 - t144 * t99) * qJD(3)) * mrSges(4,3)) * pkin(1) + t150 * t82 + t162 * t59 + (-Ifges(3,5) * t137 - Ifges(3,6) * t145) * qJD(2) + (t160 * t59 - t29 * t85 - t30 * t86 - t60 * t65) * mrSges(5,3) + (-t11 * t43 + t12 * t44 + t32 * t39 - t33 * t40) * mrSges(11,3) + t146 + t81 * t5 - t60 * t36 + t170; -0.2e1 * t188 - 0.2e1 * t204 + 0.2e1 * t31 - 0.2e1 * t48 + 0.2e1 * t53 + 0.2e1 * t54 + 0.2e1 * t56 + 0.2e1 * t72 + 0.2e1 * m(11) * (t32 * t44 + t33 * t43) + 0.2e1 * m(6) * (t168 * t82 - t60 * t81) + 0.2e1 * m(5) * (t59 * t86 + t60 * t85) + t154 + t87 + 0.2e1 * (t151 + t152) * pkin(1); t147 + ((-t135 * t30 - t143 * t29) * mrSges(5,3) + ((mrSges(5,3) * t65 + t36) * t135 + (mrSges(5,3) * t160 + t162) * t143) * qJD(4)) * pkin(5) + t150 * t120 + t121 * t5; t148 + m(6) * (-t121 * t60 + (t193 + t194) * t120) + (m(5) * (t135 * t59 + t143 * t60) + (m(6) * (t135 * t81 + t225 * t82) + m(5) * (-t135 * t85 + t143 * t86) + t166) * qJD(4)) * pkin(5) + pkin(1) * t152 + t220; 0.2e1 * t112 + 0.2e1 * t113 + 0.2e1 * t80 + 0.2e1 * t92 + 0.2e1 * (m(6) * (t120 * t225 + t121 * t135) + t166) * t199 + t154; -pkin(9) * t5 + pkin(11) * t150 + t149; m(6) * (pkin(9) * t60 + pkin(11) * t168) - t209 + t148; -t209 + (m(6) * (-pkin(9) * t135 + pkin(11) * t225) + t166) * t199 + t154 + t220; t154 - 0.2e1 * t209; t8 * t198 + (-mrSges(6,2) * t8 - Ifges(6,6) * t29) * t134 + ((-mrSges(6,2) * t34 - Ifges(6,6) * t65) * t142 + (-mrSges(6,1) * t34 - Ifges(6,5) * t65) * t134) * qJD(5) + t203; t122 - t165 * t59 + (-t82 * t198 + (mrSges(6,2) * t82 - Ifges(6,6)) * t134) * qJD(5); t122 - t165 * t176 + (t108 * t120 - t195) * qJD(5); t122 + (pkin(11) * t108 - t195) * qJD(5); 0; (Ifges(7,5) * t141 - Ifges(7,6) * t133) * qJD(6); 0; 0; 0; 0; 0; (-t11 * t57 + t12 * t58 + t39 * t50 - t40 * t51) * mrSges(11,3) + t170; m(11) * (t32 * t58 + t33 * t57 + t43 * t51 + t44 * t50) + t49 + t31 + (-t50 - t32) * mrSges(11,2) + pkin(1) * t151; 0; 0; 0; 0; 0.2e1 * m(11) * (t50 * t58 + t51 * t57) + 0.2e1 * t169; t146; t87; 0; 0; 0; 0; 0; t87; t202; t94; 0; 0; 0; 0; 0; t94; 0; t211; t31 - t188; 0; 0; 0; 0; t169; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_13_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46), t1(56), t1(67), t1(79); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47), t1(57), t1(68), t1(80); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48), t1(58), t1(69), t1(81); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49), t1(59), t1(70), t1(82); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50), t1(60), t1(71), t1(83); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51), t1(61), t1(72), t1(84); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52), t1(62), t1(73), t1(85); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53), t1(63), t1(74), t1(86); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54), t1(64), t1(75), t1(87); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55), t1(65), t1(76), t1(88); t1(56), t1(57), t1(58), t1(59), t1(60), t1(61), t1(62), t1(63), t1(64), t1(65), t1(66), t1(77), t1(89); t1(67), t1(68), t1(69), t1(70), t1(71), t1(72), t1(73), t1(74), t1(75), t1(76), t1(77), t1(78), t1(90); t1(79), t1(80), t1(81), t1(82), t1(83), t1(84), t1(85), t1(86), t1(87), t1(88), t1(89), t1(90), t1(91);];
Mq = res;
