% Calculate time derivative of joint inertia matrix for
% palh1m2OL
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh1m2OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(13,1),zeros(20,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_inertiaDJ_slag_vp2: qJ has to be [13x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [13 1]), ...
  'palh1m2OL_inertiaDJ_slag_vp2: qJD has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_inertiaDJ_slag_vp2: pkin has to be [20x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'palh1m2OL_inertiaDJ_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'palh1m2OL_inertiaDJ_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'palh1m2OL_inertiaDJ_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:17:44
% EndTime: 2020-05-02 21:18:06
% DurationCPUTime: 4.14s
% Computational Cost: add. (3242->327), mult. (7725->547), div. (0->0), fcn. (7623->20), ass. (0->175)
t138 = sin(qJ(5));
t129 = t138 ^ 2;
t146 = cos(qJ(5));
t130 = t146 ^ 2;
t191 = t129 + t130;
t187 = qJD(5) * t146;
t139 = sin(qJ(4));
t147 = cos(qJ(4));
t140 = sin(qJ(3));
t149 = cos(qJ(2));
t141 = sin(qJ(2));
t148 = cos(qJ(3));
t194 = t141 * t148;
t100 = t140 * t149 + t194;
t158 = t100 * qJD(2);
t151 = -qJD(3) * t100 - t158;
t192 = t148 * t149;
t101 = -t140 * t141 + t192;
t64 = -t100 * t139 + t101 * t147;
t236 = qJD(2) + qJD(3);
t77 = t236 * t101;
t30 = qJD(4) * t64 + t139 * t151 + t147 * t77;
t65 = t100 * t147 + t101 * t139;
t164 = t138 * t30 + t65 * t187;
t134 = sin(qJ(9));
t142 = cos(qJ(9));
t93 = (mrSges(10,1) * t134 + mrSges(10,2) * t142) * qJD(9) * pkin(2);
t86 = 0.2e1 * t93;
t123 = pkin(1) * t141 - pkin(15);
t234 = 0.2e1 * t123;
t131 = sin(pkin(19));
t132 = cos(pkin(19));
t133 = cos(qJ(10));
t213 = sin(qJ(10));
t95 = -t131 * t213 - t132 * t133;
t233 = t191 * t147;
t232 = qJD(2) + qJD(7);
t231 = qJD(2) + qJD(8);
t230 = qJD(3) + qJD(4);
t211 = pkin(5) * qJD(4);
t183 = t147 * t211;
t174 = mrSges(6,3) * t183;
t112 = t129 * t174;
t113 = t130 * t174;
t172 = mrSges(6,1) * t138 + mrSges(6,2) * t146;
t102 = t172 * qJD(5);
t126 = -pkin(5) * t147 - pkin(9);
t78 = t126 * t102;
t210 = mrSges(6,1) * t146;
t107 = mrSges(6,2) * t138 - t210;
t91 = t139 * t107 * t211;
t229 = t112 + t113 + t78 + t91;
t29 = qJD(4) * t65 + t139 * t77 - t147 * t151;
t124 = pkin(5) * t140 + pkin(1);
t189 = qJD(3) * t140;
t190 = qJD(2) * t149;
t60 = t124 * t190 + (t149 * t189 + t194 * t236) * pkin(5);
t8 = pkin(9) * t29 - pkin(11) * t30 + t60;
t228 = 0.2e1 * t8;
t227 = -0.2e1 * pkin(15);
t136 = sin(qJ(7));
t144 = cos(qJ(7));
t166 = t136 * t141 - t144 * t149;
t225 = -0.2e1 * t166;
t161 = -t131 * t133 + t132 * t213;
t98 = -t136 * t149 - t141 * t144;
t75 = t232 * t98;
t76 = t232 * t166;
t82 = t95 * qJD(10);
t83 = t161 * qJD(10);
t11 = -t161 * t76 - t166 * t83 + t75 * t95 + t82 * t98;
t12 = t161 * t75 + t166 * t82 + t76 * t95 + t83 * t98;
t222 = Ifges(11,5) * t11 + Ifges(11,6) * t12;
t221 = pkin(9) * t102;
t220 = Ifges(6,5) * t29;
t219 = Ifges(6,5) * t64;
t218 = Ifges(6,6) * t64;
t217 = t29 * Ifges(6,6);
t195 = t139 * t148;
t58 = t183 + t230 * pkin(1) * (t140 * t147 + t195);
t216 = t58 * mrSges(5,2);
t202 = t146 * t30;
t215 = -Ifges(6,5) * t202 - Ifges(6,3) * t29;
t135 = sin(qJ(8));
t143 = cos(qJ(8));
t167 = t135 * t141 - t143 * t149;
t96 = -t135 * t149 - t141 * t143;
t62 = -t134 * t167 - t142 * t96;
t73 = t231 * t96;
t74 = t231 * t167;
t27 = qJD(9) * t62 - t134 * t74 - t142 * t73;
t169 = t134 * t96 - t142 * t167;
t28 = qJD(9) * t169 + t134 * t73 - t142 * t74;
t214 = Ifges(10,5) * t27 + Ifges(10,6) * t28;
t212 = pkin(1) * qJD(7);
t209 = Ifges(6,4) * t138;
t208 = Ifges(6,4) * t146;
t207 = Ifges(6,6) * t138;
t206 = t129 * t58;
t205 = t130 * t58;
t203 = t138 * t65;
t201 = t146 * t65;
t110 = pkin(1) * t136 + pkin(4) * t131;
t111 = pkin(1) * t144 + pkin(4) * t132;
t32 = t110 * t83 + t111 * t82 + (t136 * t161 + t144 * t95) * t212;
t200 = t32 * mrSges(11,2);
t188 = qJD(5) * t138;
t79 = -pkin(5) * t192 + t124 * t141 - pkin(15);
t34 = -pkin(9) * t64 - pkin(11) * t65 + t79;
t186 = 0.2e1 * t34;
t185 = pkin(1) * t147 * t148;
t181 = pkin(1) * t140 + pkin(5);
t170 = -Ifges(6,2) * t138 + t208;
t18 = t170 * t65 - t218;
t180 = -t18 + t218;
t178 = -t188 / 0.2e1;
t177 = Ifges(8,5) * t75 + Ifges(8,6) * t76 + t222;
t51 = (-t131 * t82 + t132 * t83) * pkin(4);
t49 = t51 * mrSges(11,1);
t50 = (t131 * t83 + t132 * t82) * pkin(4);
t176 = -t50 * mrSges(11,2) + t49;
t175 = t191 * t58;
t115 = t181 * t139;
t84 = t115 - t185;
t173 = -t139 * mrSges(5,1) - t147 * mrSges(5,2);
t171 = Ifges(6,1) * t146 - t209;
t37 = mrSges(6,2) * t64 - mrSges(6,3) * t203;
t38 = -mrSges(6,1) * t64 - mrSges(6,3) * t201;
t168 = -t138 * t38 + t146 * t37;
t163 = t188 * t65 - t202;
t103 = t170 * qJD(5);
t104 = t171 * qJD(5);
t108 = Ifges(6,2) * t146 + t209;
t109 = Ifges(6,1) * t138 + t208;
t162 = t103 * t146 + t104 * t138 - t108 * t188 + t109 * t187;
t160 = qJD(3) * (mrSges(4,1) * t148 - mrSges(4,2) * t140);
t159 = qJD(7) * (-mrSges(8,1) * t136 - mrSges(8,2) * t144);
t157 = t101 * Ifges(4,2) * t100;
t156 = t123 * t100 * mrSges(4,1);
t85 = pkin(1) * t195 + t147 * t181;
t6 = mrSges(6,1) * t29 + mrSges(6,3) * t163;
t7 = -mrSges(6,2) * t29 - mrSges(6,3) * t164;
t155 = -t138 * t6 + t146 * t7 + (-t138 * t37 - t146 * t38) * qJD(5);
t127 = Ifges(6,5) * t187;
t19 = t171 * t65 - t219;
t3 = -Ifges(6,4) * t163 - Ifges(6,2) * t164 + t217;
t4 = -Ifges(6,1) * t163 - Ifges(6,4) * t164 + t220;
t154 = t138 * t4 / 0.2e1 + t18 * t178 + t19 * t187 / 0.2e1 + t29 * (Ifges(6,5) * t138 + Ifges(6,6) * t146) / 0.2e1 + t146 * t3 / 0.2e1 + Ifges(5,5) * t30 - t103 * t203 / 0.2e1 + t104 * t201 / 0.2e1 - t64 * (-Ifges(6,6) * t188 + t127) / 0.2e1 - Ifges(5,6) * t29 - t164 * t108 / 0.2e1 + (t202 / 0.2e1 + t65 * t178) * t109;
t59 = -pkin(1) * t139 * t189 - qJD(4) * t115 + t185 * t230;
t48 = t59 * t107;
t53 = mrSges(6,3) * t206;
t54 = mrSges(6,3) * t205;
t55 = t59 * mrSges(5,1);
t81 = -pkin(9) - t85;
t72 = t81 * t102;
t153 = t162 - t48 + t53 + t54 + t55 + t72 - t216;
t152 = Ifges(4,5) * t77 + Ifges(4,6) * t151 + t154;
t150 = Ifges(9,5) * t73 + Ifges(9,6) * t74 + t214 + (t142 * t27 - t134 * t28 + (t134 * t169 - t142 * t62) * qJD(9)) * pkin(2) * mrSges(10,3);
t145 = cos(qJ(6));
t137 = sin(qJ(6));
t125 = pkin(5) * t139 + pkin(11);
t87 = -pkin(2) * t96 - pkin(15);
t80 = pkin(11) + t84;
t57 = (t131 * t95 - t132 * t161) * pkin(4);
t56 = (t131 * t161 + t132 * t95) * pkin(4);
t44 = t110 * t95 - t111 * t161;
t43 = t110 * t161 + t111 * t95;
t40 = -t161 * t98 - t166 * t95;
t39 = -t161 * t166 + t95 * t98;
t36 = t172 * t65;
t35 = pkin(1) * t190 + (-t131 * t75 - t132 * t76) * pkin(4);
t33 = -t110 * t82 + t111 * t83 + (-t136 * t95 + t144 * t161) * t212;
t31 = t33 * mrSges(11,1);
t5 = mrSges(6,1) * t164 - mrSges(6,2) * t163;
t1 = [(-mrSges(9,1) * t74 + mrSges(9,2) * t73) * t227 + (t37 * t228 + (-qJD(5) * t38 + t7) * t186 + t180 * t30 + (-t217 - t3 + (-t19 + t219) * qJD(5)) * t65) * t138 + (t19 * t30 + t38 * t228 + (qJD(5) * t37 + t6) * t186 + (qJD(5) * t180 + t220 + t4) * t65) * t146 + (-0.2e1 * t60 * mrSges(5,1) + (-(2 * Ifges(5,2)) - Ifges(6,3)) * t29 + t215) * t64 + 0.2e1 * (pkin(14) * (mrSges(7,1) * t137 + mrSges(7,2) * t145) + (-t137 ^ 2 + t145 ^ 2) * Ifges(7,4) + (Ifges(7,1) - Ifges(7,2)) * t137 * t145) * qJD(6) + 0.2e1 * t62 * Ifges(10,2) * t28 + 0.2e1 * t76 * Ifges(8,2) * t98 + 0.2e1 * t39 * Ifges(11,2) * t12 + 0.2e1 * t30 * t65 * Ifges(5,1) + 0.2e1 * t11 * t40 * Ifges(11,1) + 0.2e1 * (-t157 + t156) * qJD(3) + m(6) * t191 * t34 * t228 + 0.2e1 * (m(5) * t79 + mrSges(5,2) * t65) * t60 + 0.2e1 * (mrSges(4,2) * t123 + Ifges(4,1) * t100) * t77 + 0.2e1 * (t11 * t39 + t12 * t40) * Ifges(11,4) + 0.2e1 * (-t29 * t65 + t30 * t64) * Ifges(5,4) + 0.2e1 * (t100 * t151 + t101 * t77) * Ifges(4,4) + 0.2e1 * t79 * (mrSges(5,1) * t29 + mrSges(5,2) * t30) + 0.2e1 * t87 * (-mrSges(10,1) * t28 + mrSges(10,2) * t27) + 0.2e1 * t35 * (-mrSges(11,1) * t39 + mrSges(11,2) * t40) + 0.2e1 * (m(11) * t35 - mrSges(11,1) * t12 + mrSges(11,2) * t11) * ((t131 * t166 - t132 * t98) * pkin(4) + t123) + (-0.2e1 * t157 + 0.2e1 * t156 + (mrSges(3,1) * t149 - mrSges(3,2) * t141) * t227 + 0.2e1 * (t141 ^ 2 - t149 ^ 2) * Ifges(3,4) + (-0.2e1 * mrSges(4,1) * t101 - 0.2e1 * mrSges(8,1) * t98 + 0.2e1 * mrSges(4,2) * t100 + mrSges(8,2) * t225 + (m(8) + m(4)) * t234) * t149 * pkin(1) + 0.2e1 * (Ifges(3,2) - Ifges(3,1)) * t141 * t149) * qJD(2) + 0.2e1 * (-t169 * t28 + t27 * t62) * Ifges(10,4) + 0.2e1 * ((-m(10) * t87 + mrSges(10,1) * t62 + mrSges(10,2) * t169) * pkin(2) + Ifges(9,2) * t96) * t74 + 0.2e1 * (-t166 * t76 + t75 * t98) * Ifges(8,4) + 0.2e1 * (-t167 * t74 + t73 * t96) * Ifges(9,4) - 0.2e1 * t27 * t169 * Ifges(10,1) - 0.2e1 * t73 * t167 * Ifges(9,1) + t75 * Ifges(8,1) * t225 + (-mrSges(8,1) * t76 + mrSges(8,2) * t75) * t234; t150 + t152 + t155 * t80 + ((t136 * t76 - t144 * t75 + (-t136 * t166 + t144 * t98) * qJD(7)) * mrSges(8,3) + (t101 * t189 - t140 * t77 + t148 * t158) * mrSges(4,3)) * pkin(1) + t168 * t58 + (-Ifges(3,5) * t141 - Ifges(3,6) * t149) * qJD(2) + (-t29 * t84 - t30 * t85 + t58 * t64 - t59 * t65) * mrSges(5,3) + (-t11 * t43 + t12 * t44 + t32 * t39 - t33 * t40) * mrSges(11,3) + t81 * t5 - t59 * t36 + t177; -0.2e1 * t200 - 0.2e1 * t216 + 0.2e1 * t31 - 0.2e1 * t48 + 0.2e1 * t53 + 0.2e1 * t54 + 0.2e1 * t55 + 0.2e1 * t72 + 0.2e1 * m(11) * (t32 * t44 + t33 * t43) + 0.2e1 * m(6) * (t175 * t80 - t59 * t81) + 0.2e1 * m(5) * (t58 * t84 + t59 * t85) + t162 + t86 + 0.2e1 * (t159 + t160) * pkin(1); t152 + ((-t139 * t29 - t147 * t30) * mrSges(5,3) + ((mrSges(5,3) * t65 + t36) * t139 + (mrSges(5,3) * t64 + t168) * t147) * qJD(4)) * pkin(5) + t155 * t125 + t126 * t5; t153 + m(6) * (-t126 * t59 + (t205 + t206) * t125) + (m(5) * (t139 * t58 + t147 * t59) + (m(6) * (t139 * t81 + t233 * t80) + m(5) * (-t139 * t85 + t147 * t84) + t173) * qJD(4)) * pkin(5) + pkin(1) * t160 + t229; 0.2e1 * t112 + 0.2e1 * t113 + 0.2e1 * t78 + 0.2e1 * t91 + 0.2e1 * (m(6) * (t125 * t233 + t126 * t139) + t173) * t211 + t162; -pkin(9) * t5 + pkin(11) * t155 + t154; m(6) * (pkin(9) * t59 + pkin(11) * t175) - t221 + t153; -t221 + (m(6) * (-pkin(9) * t139 + pkin(11) * t233) + t173) * t211 + t162 + t229; t162 - 0.2e1 * t221; t8 * t210 + (-mrSges(6,2) * t8 - Ifges(6,6) * t30) * t138 + ((-mrSges(6,2) * t34 - Ifges(6,6) * t65) * t146 + (-mrSges(6,1) * t34 - Ifges(6,5) * t65) * t138) * qJD(5) - t215; t127 - t172 * t58 + (-t80 * t210 + (mrSges(6,2) * t80 - Ifges(6,6)) * t138) * qJD(5); t127 - t172 * t183 + (t107 * t125 - t207) * qJD(5); t127 + (pkin(11) * t107 - t207) * qJD(5); 0; (Ifges(7,5) * t145 - Ifges(7,6) * t137) * qJD(6); 0; 0; 0; 0; 0; (-t11 * t56 + t12 * t57 + t39 * t50 - t40 * t51) * mrSges(11,3) + t177; t31 + t49 + m(11) * (t32 * t57 + t33 * t56 + t43 * t51 + t44 * t50) + (-t32 - t50) * mrSges(11,2) + pkin(1) * t159; 0; 0; 0; 0; 0.2e1 * m(11) * (t50 * t57 + t51 * t56) + 0.2e1 * t176; t150; t86; 0; 0; 0; 0; 0; t86; t214; t93; 0; 0; 0; 0; 0; t93; 0; t222; t31 - t200; 0; 0; 0; 0; t176; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_13_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16), t1(22), t1(29), t1(37), t1(46), t1(56), t1(67), t1(79); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17), t1(23), t1(30), t1(38), t1(47), t1(57), t1(68), t1(80); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18), t1(24), t1(31), t1(39), t1(48), t1(58), t1(69), t1(81); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19), t1(25), t1(32), t1(40), t1(49), t1(59), t1(70), t1(82); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20), t1(26), t1(33), t1(41), t1(50), t1(60), t1(71), t1(83); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21), t1(27), t1(34), t1(42), t1(51), t1(61), t1(72), t1(84); t1(22), t1(23), t1(24), t1(25), t1(26), t1(27), t1(28), t1(35), t1(43), t1(52), t1(62), t1(73), t1(85); t1(29), t1(30), t1(31), t1(32), t1(33), t1(34), t1(35), t1(36), t1(44), t1(53), t1(63), t1(74), t1(86); t1(37), t1(38), t1(39), t1(40), t1(41), t1(42), t1(43), t1(44), t1(45), t1(54), t1(64), t1(75), t1(87); t1(46), t1(47), t1(48), t1(49), t1(50), t1(51), t1(52), t1(53), t1(54), t1(55), t1(65), t1(76), t1(88); t1(56), t1(57), t1(58), t1(59), t1(60), t1(61), t1(62), t1(63), t1(64), t1(65), t1(66), t1(77), t1(89); t1(67), t1(68), t1(69), t1(70), t1(71), t1(72), t1(73), t1(74), t1(75), t1(76), t1(77), t1(78), t1(90); t1(79), t1(80), t1(81), t1(82), t1(83), t1(84), t1(85), t1(86), t1(87), t1(88), t1(89), t1(90), t1(91);];
Mq = res;
