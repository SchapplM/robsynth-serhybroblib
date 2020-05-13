% Calculate vector of inverse dynamics joint torques for
% picker2Dm1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% qJDD [12x1]
%   Generalized joint accelerations
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
% Ifges [11x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [12x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-11 05:46
% Revision: 52c9de996ed1e2eb3528e92ec0df589a9be0640a (2020-05-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = picker2Dm1OL_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),zeros(12,1),zeros(3,1),zeros(8,1),zeros(11,1),zeros(11,3),zeros(11,6)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: qJD has to be [12x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [12 1]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: qJDD has to be [12x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [11 1]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: m has to be [11x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [11,3]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: mrSges has to be [11x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [11 6]), ...
  'picker2Dm1OL_invdynJ_fixb_slag_vp2: Ifges has to be [11x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-11 05:44:58
% EndTime: 2020-05-11 05:45:08
% DurationCPUTime: 3.07s
% Computational Cost: add. (2742->382), mult. (5000->506), div. (0->0), fcn. (2594->34), ass. (0->204)
t178 = cos(qJ(4));
t131 = pkin(3) * t178 + pkin(4);
t164 = sin(qJ(10));
t165 = cos(qJ(10));
t170 = sin(qJ(4));
t206 = qJD(10) * t170;
t208 = qJD(10) * t164;
t226 = t165 * t170;
t180 = cos(qJ(2));
t172 = sin(qJ(2));
t221 = t172 * t178;
t188 = -t170 * t180 - t221;
t235 = pkin(1) * qJD(1);
t62 = t188 * t235;
t224 = t170 * t172;
t187 = t178 * t180 - t224;
t65 = t187 * t235;
t255 = -t164 * t65 + t165 * t62 + t131 * t208 + (t165 * t206 + (t164 * t178 + t226) * qJD(4)) * pkin(3);
t207 = qJD(10) * t165;
t227 = t164 * t170;
t254 = -t164 * t62 - t165 * t65 + t131 * t207 - (t164 * t206 + (-t165 * t178 + t227) * qJD(4)) * pkin(3);
t179 = cos(qJ(3));
t132 = -pkin(2) * t179 + pkin(6);
t166 = sin(qJ(9));
t171 = sin(qJ(3));
t174 = cos(qJ(9));
t210 = qJD(9) * t174;
t211 = qJD(9) * t166;
t222 = t171 * t174;
t220 = t172 * t179;
t186 = t171 * t180 + t220;
t63 = t186 * t235;
t137 = t180 * t235;
t203 = t172 * t235;
t98 = t171 * t203;
t66 = -t137 * t179 + t98;
t253 = -t166 * t66 + t174 * t63 + t132 * t211 + (-t171 * t210 + (-t166 * t179 - t222) * qJD(3)) * pkin(2);
t225 = t166 * t171;
t252 = -t166 * t63 - t174 * t66 + t132 * t210 - (-t171 * t211 + (t174 * t179 - t225) * qJD(3)) * pkin(2);
t158 = qJDD(1) + qJDD(2);
t142 = qJDD(3) + t158;
t197 = t179 * t203;
t216 = qJD(3) * t171;
t156 = t180 * pkin(1);
t86 = -qJD(2) * t203 + qJDD(1) * t156;
t58 = pkin(2) * t158 + t86;
t217 = qJD(2) * t180;
t87 = (qJD(1) * t217 + qJDD(1) * t172) * pkin(1);
t160 = qJD(1) + qJD(2);
t89 = pkin(2) * t160 + t137;
t19 = qJD(3) * t197 + t171 * t87 - t179 * t58 + t89 * t216;
t251 = t19 * mrSges(4,1) + Ifges(4,3) * t142;
t163 = qJ(1) + qJ(2);
t153 = qJ(4) + t163;
t130 = qJ(10) + t153;
t106 = sin(t130);
t107 = cos(t130);
t154 = qJ(3) + t163;
t126 = sin(t154);
t129 = cos(t154);
t149 = sin(t163);
t151 = cos(t163);
t125 = sin(t153);
t128 = cos(t153);
t219 = t125 * mrSges(5,1) + t128 * mrSges(5,2);
t136 = qJ(9) + t154;
t112 = sin(t136);
t113 = cos(t136);
t236 = t112 * mrSges(10,1) + t113 * mrSges(10,2);
t152 = qJ(6) + t163;
t124 = sin(t152);
t127 = cos(t152);
t77 = -mrSges(7,1) * t124 - mrSges(7,2) * t127;
t250 = t106 * mrSges(11,1) - t149 * mrSges(3,1) + t126 * mrSges(4,1) + t107 * mrSges(11,2) - t151 * mrSges(3,2) + t129 * mrSges(4,2) - t219 - t236 - t77;
t201 = -t113 * mrSges(10,1) + t112 * mrSges(10,2);
t195 = -t129 * mrSges(4,1) - t201;
t91 = t107 * mrSges(11,1);
t237 = -t125 * mrSges(5,2) - t91;
t243 = m(11) + m(5);
t244 = m(4) + m(10);
t78 = t127 * mrSges(7,1) - t124 * mrSges(7,2);
t249 = -t78 + t106 * mrSges(11,2) + t126 * mrSges(4,2) + (t244 * pkin(2) + pkin(3) * t243 + mrSges(3,1)) * t151 + (m(11) * pkin(4) + mrSges(5,1)) * t128 - t149 * mrSges(3,2) - m(10) * pkin(6) * t129 + t195 + t237;
t88 = pkin(3) * t160 + t137;
t52 = -t170 * t203 + t178 * t88;
t143 = qJD(6) + t160;
t169 = sin(qJ(6));
t177 = cos(qJ(6));
t189 = t169 * t172 - t177 * t180;
t82 = t189 * pkin(1);
t61 = qJD(1) * t82;
t29 = qJD(6) * t61 - t169 * t86 - t177 * t87;
t190 = t169 * t180 + t172 * t177;
t64 = t190 * t235;
t248 = -t64 * t143 * mrSges(7,1) + (t143 * t61 - t29) * mrSges(7,2);
t247 = pkin(1) * (qJD(2) + qJD(6));
t246 = pkin(2) * m(4);
t245 = pkin(3) * m(5);
t141 = qJDD(4) + t158;
t108 = qJDD(10) + t141;
t54 = t170 * t88 + t178 * t203;
t57 = pkin(3) * t158 + t86;
t17 = -t54 * qJD(4) - t170 * t87 + t178 * t57;
t12 = pkin(4) * t141 + t17;
t16 = t52 * qJD(4) + t170 * t57 + t178 * t87;
t231 = t165 * t54;
t144 = qJD(4) + t160;
t49 = pkin(4) * t144 + t52;
t21 = t164 * t49 + t231;
t4 = qJD(10) * t21 - t12 * t165 + t16 * t164;
t242 = t4 * mrSges(11,1) + Ifges(11,3) * t108;
t114 = qJDD(9) + t142;
t13 = pkin(6) * t142 + t19;
t53 = -t179 * t89 + t98;
t18 = t53 * qJD(3) - t171 * t58 - t179 * t87;
t55 = t171 * t89 + t197;
t229 = t174 * t55;
t145 = qJD(3) + t160;
t50 = pkin(6) * t145 + t53;
t23 = t166 * t50 - t229;
t6 = qJD(9) * t23 - t13 * t174 + t166 * t18;
t241 = t6 * mrSges(10,1) + Ifges(10,3) * t114;
t240 = g(1) * t106;
t239 = g(1) * t126;
t117 = qJD(9) + t145;
t234 = mrSges(10,1) * t117;
t233 = t144 * mrSges(5,1);
t232 = t164 * t54;
t230 = t166 * t55;
t140 = qJDD(6) + t158;
t30 = qJD(6) * t64 + t169 * t87 - t177 * t86;
t228 = t30 * mrSges(7,1) + Ifges(7,3) * t140;
t223 = t171 * t172;
t122 = pkin(3) * t149;
t173 = sin(qJ(1));
t155 = t173 * pkin(1);
t218 = t122 + t155;
t215 = qJD(3) * t179;
t214 = qJD(4) * t170;
t213 = qJD(4) * t178;
t159 = qJD(1) + qJD(8);
t212 = qJD(8) * t159;
t209 = -m(3) - m(7) - m(9);
t205 = qJD(1) * qJD(8);
t123 = pkin(2) * t149;
t202 = -pkin(6) * t126 + t123;
t200 = t17 * mrSges(5,1) + Ifges(5,3) * t141 + t242;
t162 = qJ(1) + qJ(8);
t148 = sin(t162);
t150 = cos(t162);
t199 = t150 * mrSges(9,1) - t148 * mrSges(9,2);
t134 = t156 + pkin(2);
t72 = pkin(1) * t223 - t134 * t179;
t133 = t156 + pkin(3);
t71 = -pkin(1) * t224 + t178 * t133;
t157 = qJDD(1) + qJDD(8);
t167 = sin(qJ(8));
t175 = cos(qJ(8));
t84 = (-qJDD(1) * t175 + t167 * t205) * pkin(1);
t85 = (-qJDD(1) * t167 - t175 * t205) * pkin(1);
t196 = t84 * mrSges(9,1) - t85 * mrSges(9,2) + Ifges(9,3) * t157;
t194 = -mrSges(9,1) * t148 - mrSges(9,2) * t150;
t20 = -t165 * t49 + t232;
t59 = pkin(4) + t71;
t74 = pkin(1) * t221 + t133 * t170;
t33 = t164 * t74 - t165 * t59;
t193 = t164 * t59 + t165 * t74;
t22 = -t174 * t50 - t230;
t60 = pkin(6) + t72;
t75 = -pkin(1) * t220 - t134 * t171;
t37 = t166 * t75 - t174 * t60;
t192 = t166 * t60 + t174 * t75;
t191 = -g(1) * t236 + t241;
t3 = qJD(10) * t20 - t12 * t164 - t16 * t165;
t185 = g(1) * t107 + g(2) * t106 - t3;
t184 = t86 * mrSges(3,1) - t87 * mrSges(3,2) + Ifges(3,3) * t158 + t200 + t228 + t241 + t251;
t181 = cos(qJ(1));
t176 = cos(qJ(7));
t168 = sin(qJ(7));
t161 = pkin(8) + qJ(5);
t147 = cos(t161);
t146 = sin(t161);
t116 = qJD(10) + t144;
t105 = pkin(4) * t125;
t83 = t190 * pkin(1);
t73 = pkin(2) * t222 - t132 * t166;
t70 = -pkin(2) * t225 - t132 * t174;
t68 = -pkin(3) * t226 - t131 * t164;
t67 = pkin(3) * t227 - t131 * t165;
t48 = t190 * t247;
t47 = t189 * t247;
t46 = t134 * t216 + (qJD(2) * t186 + t172 * t215) * pkin(1);
t45 = -t134 * t215 + (t172 * t216 + (-t179 * t180 + t223) * qJD(2)) * pkin(1);
t44 = -t133 * t214 + (qJD(2) * t188 - t172 * t213) * pkin(1);
t43 = t133 * t213 + (qJD(2) * t187 - t172 * t214) * pkin(1);
t28 = -t174 * t53 - t230;
t27 = t166 * t53 - t229;
t26 = -t165 * t52 + t232;
t25 = t164 * t52 + t231;
t10 = qJD(9) * t192 + t166 * t45 - t174 * t46;
t9 = qJD(9) * t37 - t166 * t46 - t174 * t45;
t8 = qJD(10) * t193 + t164 * t43 - t165 * t44;
t7 = qJD(10) * t33 - t164 * t44 - t165 * t43;
t5 = qJD(9) * t22 - t13 * t166 - t174 * t18;
t1 = [t196 + (-m(4) * (t123 + t155) - m(5) * t218 - m(11) * (t105 + t218) - m(10) * (t155 + t202) - mrSges(2,2) * t181 - t173 * mrSges(2,1) - t194 + t250) * g(1) + (mrSges(2,1) * t181 - t173 * mrSges(2,2) - t199 + t249) * g(2) + t184 + m(5) * (t16 * t74 + t17 * t71 + t54 * t43 + t52 * t44) + m(4) * (t18 * t75 + t19 * t72 - t55 * t45 + t53 * t46) + m(7) * (-t29 * t83 + t30 * t82 - t47 * t64 + t48 * t61) + m(11) * (-t193 * t3 + t20 * t8 - t21 * t7 + t33 * t4) + m(10) * (t10 * t22 - t192 * t5 - t23 * t9 + t37 * t6) + (t82 * t140 + t48 * t143) * mrSges(7,1) + (t83 * t140 - t47 * t143 - t29) * mrSges(7,2) + (-t75 * t142 - t45 * t145 - t18) * mrSges(4,2) + (t108 * t193 - t7 * t116 - t3) * mrSges(11,2) + (t10 * t117 + t37 * t114) * mrSges(10,1) + (t72 * t142 + t46 * t145) * mrSges(4,1) + (t71 * t141 + t44 * t144) * mrSges(5,1) + (t33 * t108 + t8 * t116) * mrSges(11,1) + (t114 * t192 - t9 * t117 - t5) * mrSges(10,2) + (-t74 * t141 - t43 * t144 - t16) * mrSges(5,2) + Ifges(2,3) * qJDD(1) + (t209 * t173 * g(1) + (-t209 + t243 + t244) * t181 * g(2) + (t157 * t167 + t175 * t212) * mrSges(9,2) + (-t158 * t172 - t160 * t217) * mrSges(3,2) + (-t157 * t175 + t167 * t212) * mrSges(9,1) + (-qJD(2) * t160 * t172 + t158 * t180) * mrSges(3,1) + m(9) * (-t167 * t85 - t175 * t84) + m(3) * (t172 * t87 + t180 * t86)) * pkin(1); (t3 * t68 + t4 * t67 + (-t122 - t105) * g(1) + t254 * t21 + t255 * t20) * m(11) + (-t202 * g(1) + t253 * t22 + t252 * t23 + t5 * t73 + t6 * t70) * m(10) - t62 * t233 + ((t142 * t171 + t145 * t215) * mrSges(4,2) + (-t142 * t179 + t145 * t216) * mrSges(4,1)) * pkin(2) + ((-t141 * t170 - t144 * t213) * mrSges(5,2) + (t141 * t178 - t144 * t214) * mrSges(5,1)) * pkin(3) + (mrSges(3,1) * t172 + mrSges(3,2) * t180) * t160 * t235 + (t16 * t170 + t17 * t178 + (-t170 * t52 + t178 * t54) * qJD(4)) * t245 + (-t171 * t18 - t179 * t19 + (t171 * t53 + t179 * t55) * qJD(3)) * t246 + t249 * g(2) + t184 + ((-t245 - t246) * t149 + t250) * g(1) + (-t68 * t108 + t254 * t116 - t3) * mrSges(11,2) + (-t73 * t114 + t252 * t117 - t5) * mrSges(10,2) + (t70 * t114 + t253 * t117) * mrSges(10,1) + (t67 * t108 + t255 * t116) * mrSges(11,1) - m(5) * (t52 * t62 + t54 * t65) - m(4) * (t53 * t63 - t55 * t66) - t63 * t145 * mrSges(4,1) + (t144 * t65 - t16) * mrSges(5,2) + (t145 * t66 - t18) * mrSges(4,2) + t248; -m(10) * (t22 * t27 - t23 * t28) - t27 * t234 + t195 * g(2) + (t117 * t28 - t5) * mrSges(10,2) + (-t145 * t55 + t239) * mrSges(4,1) + (g(1) * t129 + g(2) * t126 + t145 * t53 - t18) * mrSges(4,2) + ((t114 * t166 + t117 * t210) * mrSges(10,2) + (-t114 * t174 + t117 * t211) * mrSges(10,1) + (-g(2) * t129 - t166 * t5 - t174 * t6 + t210 * t23 + t211 * t22 + t239) * m(10)) * pkin(6) + t191 + t251; -m(11) * (t20 * t25 - t21 * t26) + t54 * t233 - g(1) * t219 + (t128 * mrSges(5,1) + t237) * g(2) + (t144 * t52 - t16) * mrSges(5,2) + (-t116 * t25 + t240) * mrSges(11,1) + (t116 * t26 + t185) * mrSges(11,2) + ((t108 * t164 + t116 * t207) * mrSges(11,2) + (-t108 * t165 + t116 * t208) * mrSges(11,1) + (-g(1) * t125 + g(2) * t128 - t164 * t3 - t165 * t4 + t20 * t208 + t207 * t21) * m(11)) * pkin(4) + t200; Ifges(6,3) * qJDD(5) - g(1) * (-mrSges(6,1) * t146 - mrSges(6,2) * t147) - g(2) * (mrSges(6,1) * t147 - mrSges(6,2) * t146); -g(1) * t77 - g(2) * t78 + t228 + t248; Ifges(8,3) * qJDD(7) - g(1) * (mrSges(8,1) * t176 - mrSges(8,2) * t168) - g(2) * (mrSges(8,1) * t168 + mrSges(8,2) * t176); -g(1) * t194 - g(2) * t199 + (-mrSges(9,1) * t167 - mrSges(9,2) * t175) * t159 * t235 + t196; -t23 * t234 - g(2) * t201 + (t117 * t22 - t5) * mrSges(10,2) + t191; -g(2) * t91 + (-t116 * t21 + t240) * mrSges(11,1) + (t116 * t20 + t185) * mrSges(11,2) + t242; 0; 0;];
tau = t1;
