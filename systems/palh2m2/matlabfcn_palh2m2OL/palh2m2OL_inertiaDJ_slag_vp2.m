% Calculate time derivative of joint inertia matrix for
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2OL_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2OL_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'palh2m2OL_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'palh2m2OL_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:48:15
% EndTime: 2020-05-03 01:48:44
% DurationCPUTime: 5.28s
% Computational Cost: add. (4021->333), mult. (6605->541), div. (0->0), fcn. (4999->10), ass. (0->195)
t135 = qJD(2) + qJD(3);
t149 = sin(qJ(4));
t154 = cos(qJ(4));
t134 = qJD(3) + qJD(4);
t129 = qJD(2) + t134;
t231 = mrSges(7,1) * qJD(6);
t131 = pkin(5) * t231;
t152 = cos(qJ(6));
t136 = t152 ^ 2;
t147 = sin(qJ(6));
t148 = sin(qJ(5));
t153 = cos(qJ(5));
t144 = Ifges(7,2) - Ifges(7,1);
t208 = qJD(6) * t144;
t133 = qJD(4) + qJD(5);
t128 = qJD(3) + t133;
t121 = qJD(2) + t128;
t222 = t148 * t121;
t162 = Ifges(7,4) * t222 + t153 * t208;
t229 = mrSges(7,2) * qJD(6);
t198 = pkin(5) * t229;
t232 = mrSges(7,1) * pkin(3);
t241 = Ifges(7,6) * t121;
t246 = Ifges(7,4) * t147;
t273 = (t232 + 0.4e1 * t246) * qJD(6) + t241;
t295 = t153 * t273;
t201 = pkin(3) * qJD(6);
t102 = mrSges(7,2) * t201 - Ifges(7,5) * t121;
t277 = t102 * t147 + t208 + (t133 + t135) * Ifges(6,6);
t173 = pkin(3) * t121;
t228 = Ifges(7,6) * qJD(6);
t103 = mrSges(7,1) * t173 + t228;
t293 = t121 * (Ifges(7,4) + Ifges(6,5)) + t103 * t147;
t38 = t293 * t148 + t153 * t277;
t141 = Ifges(7,5) * qJD(6);
t101 = mrSges(7,2) * t173 - t141;
t223 = t144 * t147;
t60 = t121 * t223 + t101;
t17 = -Ifges(5,6) * t129 - t147 * t198 - (t60 * t148 - t131 - t295) * t152 + 0.2e1 * t162 * t136 - t38;
t254 = pkin(5) * t129;
t117 = mrSges(7,2) * t254;
t217 = t153 * t121;
t163 = Ifges(7,4) * t217 - t148 * t208;
t239 = t147 * mrSges(7,1);
t297 = t148 * t273;
t39 = -t148 * t277 + t293 * t153;
t271 = -0.2e1 * t163 * t136 + (t60 * t153 + t117 + t297) * t152 + t129 * (Ifges(5,5) + (-mrSges(6,3) + t239) * pkin(5)) + t39;
t6 = t149 * t271 - t17 * t154;
t305 = -Ifges(4,6) * t135 - t6;
t150 = sin(qJ(3));
t111 = t152 * mrSges(7,2) + t239;
t171 = -mrSges(5,3) - mrSges(6,3) + t111;
t7 = t17 * t149 + t154 * t271;
t2 = t135 * (t171 * pkin(2) + Ifges(4,5)) + t7;
t304 = t150 * t2;
t274 = t152 * mrSges(7,1) - t147 * mrSges(7,2);
t100 = -pkin(3) * m(7) - mrSges(6,1) - t274;
t155 = cos(qJ(3));
t227 = t128 * t150;
t146 = mrSges(6,2) + mrSges(7,3);
t215 = t147 * t231 + t152 * t229;
t84 = t128 * t146 + t215;
t184 = (t100 * t227 - t155 * t84) * pkin(4);
t87 = t146 * t133 + t215;
t216 = t87 * pkin(2);
t303 = -t184 + t216;
t302 = t100 * t148;
t202 = pkin(2) * qJD(6);
t104 = mrSges(7,1) * t201 + t241;
t266 = (t101 * t148 - t104 * t153 - t131) * t147 - (t102 * t153 + t103 * t148 + t198) * t152 - Ifges(7,3) * t222;
t29 = (t101 * t153 + t104 * t148 + t117) * t147 - (mrSges(7,1) * t254 - t102 * t148 + t103 * t153) * t152 - Ifges(7,3) * t217;
t299 = -t111 * t202 + t29 * t149 + t154 * t266;
t142 = pkin(5) * qJD(5);
t197 = t142 * t100;
t259 = 0.2e1 * t148;
t252 = pkin(5) * t153;
t97 = qJD(5) * t146 + t215;
t80 = -0.2e1 * t97 * t252;
t298 = t197 * t259 + t80;
t230 = mrSges(7,2) * pkin(3);
t179 = t223 + t230;
t169 = t179 * t152;
t193 = t147 * t232;
t247 = Ifges(7,4) * t136;
t292 = t169 + Ifges(7,4) + t193 - 0.2e1 * t247;
t226 = t128 * t155;
t113 = pkin(4) * t226;
t214 = pkin(2) * t133;
t255 = pkin(4) * t150;
t69 = t84 * t255;
t268 = t69 + (t113 + t214) * t100;
t287 = t268 * t148 - t153 * t303;
t256 = pkin(4) * t134;
t116 = t155 * t256;
t203 = pkin(2) * qJD(4);
t108 = t116 + t203;
t114 = mrSges(5,2) * t134 * t255;
t258 = m(6) + m(7);
t183 = t258 * pkin(5) + mrSges(5,1);
t286 = (-t108 * t183 + t148 * t303 + t153 * t268 + t114) * t149;
t285 = 0.2e1 * pkin(4);
t151 = sin(qJ(2));
t156 = cos(qJ(2));
t105 = t150 * t156 + t155 * t151;
t219 = t150 * t151;
t106 = t155 * t156 - t219;
t62 = t154 * t105 + t106 * t149;
t63 = t135 * t105;
t64 = t135 * t106;
t27 = -t62 * qJD(4) - t149 * t64 - t154 * t63;
t61 = -t149 * t105 + t154 * t106;
t28 = t61 * qJD(4) - t63 * t149 + t154 * t64;
t34 = t61 * t148 + t153 * t62;
t8 = -t34 * qJD(5) - t148 * t28 + t153 * t27;
t284 = 0.2e1 * t8;
t282 = -0.2e1 * qJD(6);
t281 = 0.2e1 * qJD(6);
t74 = t149 * t216;
t267 = t74 + (t154 * t214 + t142) * t100;
t122 = t154 * pkin(5) + pkin(2);
t210 = qJD(3) * t155;
t211 = qJD(3) * t150;
t212 = qJD(2) * t156;
t213 = qJD(2) * t151;
t218 = t150 * t154;
t220 = t150 * t149;
t221 = t149 * t155;
t95 = pkin(5) * t220 - t122 * t155 - pkin(4);
t98 = pkin(5) * t221 + t150 * t122;
t23 = -t95 * t213 + (t122 * t211 + (t149 * t210 + (t218 + t221) * qJD(4)) * pkin(5)) * t156 + (t122 * t210 + (-t149 * t211 + (t154 * t155 - t220) * qJD(4)) * pkin(5)) * t151 + t98 * t212;
t5 = -t8 * pkin(3) + t23;
t265 = -0.2e1 * t5;
t264 = 2 * Ifges(6,4);
t261 = -0.2e1 * t156 * pkin(4) - (2 * pkin(1));
t260 = -0.2e1 * t147;
t257 = pkin(5) * t97;
t253 = pkin(5) * t148;
t251 = mrSges(5,2) * t154;
t250 = mrSges(7,3) * t147;
t249 = mrSges(7,3) * t152;
t248 = Ifges(7,1) * t152;
t245 = Ifges(7,4) * t152;
t243 = Ifges(7,2) * t147;
t225 = t133 * t149;
t176 = t100 * t225 - t154 * t87;
t234 = (t176 * pkin(2) - t257) * t153;
t224 = t133 * t154;
t209 = qJD(6) * t111;
t207 = qJD(6) * t149;
t206 = qJD(6) * t150;
t205 = qJD(6) * t154;
t204 = qJD(6) * t155;
t200 = 0.2e1 * t156;
t66 = t292 * t282;
t199 = t66 + t298;
t196 = t153 * t142;
t48 = t98 * t151 + t95 * t156 - pkin(1);
t187 = -0.2e1 * pkin(3) - t252;
t186 = mrSges(7,2) * t196;
t185 = t183 * t134 * t218;
t181 = Ifges(7,5) * t152 - Ifges(7,6) * t147;
t33 = -t148 * t62 + t153 * t61;
t180 = t149 * t162 - t154 * t163;
t172 = t274 * qJD(6);
t168 = (-mrSges(4,2) * t155 - ((m(5) + t258) * pkin(2) + mrSges(4,1)) * t150) * qJD(3);
t167 = (-t149 * t183 - t251) * qJD(4);
t132 = mrSges(7,1) * t142;
t123 = -t155 * pkin(2) - pkin(4);
t107 = t116 + 0.2e1 * t203;
t79 = -mrSges(7,1) * t202 + (-mrSges(7,1) * t204 + mrSges(7,2) * t227) * pkin(4);
t78 = -mrSges(7,2) * t202 + (-mrSges(7,1) * t227 - mrSges(7,2) * t204) * pkin(4);
t71 = mrSges(7,1) * t214 + (mrSges(7,1) * t226 - mrSges(7,2) * t206) * pkin(4);
t70 = mrSges(7,2) * t214 + (mrSges(7,1) * t206 + mrSges(7,2) * t226) * pkin(4);
t68 = t179 * t121 - t141;
t65 = t292 * t281;
t56 = -t123 * t213 + (t105 * qJD(3) + t150 * t212) * pkin(2);
t47 = t149 * t163 + t154 * t162;
t46 = -t68 * t148 + t295;
t45 = t68 * t153 + t297;
t43 = t69 + t100 * (t113 + 0.2e1 * t214);
t35 = t184 - 0.2e1 * t216;
t26 = -t45 * t149 + t46 * t154;
t25 = t149 * t46 + t45 * t154;
t22 = -t149 * t38 + t39 * t154;
t21 = t149 * t39 + t38 * t154;
t20 = t33 * mrSges(7,1) - t34 * t249;
t19 = -t33 * mrSges(7,2) - t34 * t250;
t12 = -pkin(2) * t135 * t274 - t149 * t266 + t29 * t154;
t11 = Ifges(7,5) * t33 + (-t246 + t248) * t34;
t10 = Ifges(7,6) * t33 + (-t243 + t245) * t34;
t9 = t33 * qJD(5) + t27 * t148 + t153 * t28;
t4 = -t202 * t274 - t305;
t3 = pkin(2) * t172 + t305;
t1 = t2 * t155;
t13 = [(t63 * mrSges(4,1) + t64 * mrSges(4,2)) * t261 + 0.2e1 * t105 * t64 * Ifges(4,1) - 0.2e1 * t63 * Ifges(4,2) * t106 + 0.2e1 * t61 * Ifges(5,2) * t27 + 0.2e1 * t62 * t28 * Ifges(5,1) + 0.2e1 * t56 * (-t61 * mrSges(5,1) + t62 * mrSges(5,2)) + (t9 * t11 + t20 * t265) * t152 + (-t9 * t10 + t19 * t265) * t147 + (-0.2e1 * t23 * mrSges(6,1) + (t264 + t181) * t9) * t33 + ((-(pkin(1) * mrSges(3,2)) + Ifges(3,4) * t156) * t200 + (m(4) * pkin(4) * t261 + (-t106 * mrSges(4,1) + t105 * mrSges(4,2)) * t285 - (2 * pkin(1) * mrSges(3,1)) - 0.2e1 * Ifges(3,4) * t151 + (Ifges(3,1) - Ifges(3,2)) * t200) * t151) * qJD(2) + (0.2e1 * t23 * mrSges(6,2) + 0.2e1 * Ifges(6,1) * t9 + t8 * t264 + (Ifges(7,5) * t284 + t9 * t248) * t152 + (-0.2e1 * Ifges(7,6) * t8 + (t243 - 0.2e1 * t245) * t9) * t147 + (-t147 * t11 - t152 * t10 + t33 * (-Ifges(7,5) * t147 - Ifges(7,6) * t152) + (t152 * (-Ifges(7,1) * t147 - t245) - t147 * (-Ifges(7,2) * t152 - t246)) * t34) * qJD(6)) * t34 + 0.2e1 * (t62 * t27 + t61 * t28) * Ifges(5,4) + 0.2e1 * (-t105 * t63 + t64 * t106) * Ifges(4,4) + (Ifges(6,2) + Ifges(7,3)) * t33 * t284 + 0.2e1 * (m(6) * t23 - t8 * mrSges(6,1) + t9 * mrSges(6,2)) * t48 + 0.2e1 * (m(5) * t56 - t27 * mrSges(5,1) + t28 * mrSges(5,2)) * (pkin(2) * t219 + t123 * t156 - pkin(1)) + ((-t8 * mrSges(7,2) - t9 * t250) * t260 - 0.2e1 * t152 * (t8 * mrSges(7,1) - t9 * t249) + 0.2e1 * m(7) * (t147 ^ 2 + t136) * t5 + (t147 * t20 - t152 * t19) * t281) * (-t33 * pkin(3) + t48); (qJD(2) * Ifges(3,5) + t3 * t150 + t1) * t156 + t151 * (-qJD(2) * Ifges(3,6) + t3 * t155 - t304) + ((-mrSges(4,3) + t171) * t212 + t151 * t172) * pkin(4); -0.2e1 * t185 * pkin(4) - 0.2e1 * t108 * t251 + 0.2e1 * t287 * t154 + t168 * t285 + t199 + 0.2e1 * t286; (-t150 * t4 + t1) * t156 - t151 * (t4 * t155 + t304); (-t107 * t183 - t35 * t148 + t43 * t153 + t114) * t149 + (-mrSges(5,2) * t107 + t43 * t148 + t35 * t153) * t154 + (t168 - t185) * pkin(4) + t199; 0.2e1 * pkin(2) * t167 + t267 * t259 + 0.2e1 * t234 - t65; (-t6 * t150 + t7 * t155) * t156 - t151 * (t150 * t7 + t6 * t155); t286 + (-mrSges(5,2) * t203 + (-mrSges(5,2) * t155 - t150 * t183) * t256 + t287) * t154 + t199; t80 + (t74 + 0.2e1 * t197) * t148 - t65 + (t153 * t176 + t224 * t302 + t167) * pkin(2); (0.2e1 * t247 + (t232 + t246) * t260) * qJD(6) + t169 * t282 + t298; 0.2e1 * ((t47 * t150 + t155 * t180) * t156 + (-t150 * t180 + t47 * t155) * t151) * t136 + ((t150 * t26 + t25 * t155) * t156 + (-t25 * t150 + t26 * t155) * t151) * t152 + (-t150 * t21 + t22 * t155) * t156 - (t150 * t22 + t21 * t155) * t151; (t149 * t268 - t154 * t303 - t257) * t153 + (t149 * t303 + t154 * t268 + t197) * t148 + 0.4e1 * (((-Ifges(7,2) / 0.2e1 + Ifges(7,1) / 0.2e1) * t147 - t230 / 0.2e1) * t152 - t193 / 0.2e1 + (t136 - 0.1e1 / 0.2e1) * Ifges(7,4)) * qJD(6); t267 * t148 + t234 + t66; (-t146 * t153 + t302) * t142 + ((mrSges(7,2) * t187 - 0.2e1 * t223) * t152 + (0.4e1 * t136 - 0.2e1) * Ifges(7,4) + t187 * t239) * qJD(6); t66; (-pkin(4) * t209 + t12 * t150 + t155 * t299) * t156 + (-pkin(4) * qJD(2) * t274 + t12 * t155 - t150 * t299) * t151 - pkin(1) * t209; ((-mrSges(7,2) * t142 + t79 * t149 - t70 * t154) * t153 + (t70 * t149 + t79 * t154 - t131) * t148 - t141) * t152 - ((t78 * t149 + t71 * t154 + t132) * t153 + (-t71 * t149 + t78 * t154 - t198) * t148 - t228) * t147; (-t131 * t148 - t141 - t186) * t152 - t147 * (t132 * t153 - t148 * t198 - t228) + (((-mrSges(7,1) * t207 - mrSges(7,2) * t224) * t153 + (-mrSges(7,1) * t205 + mrSges(7,2) * t225) * t148) * t152 - t147 * ((mrSges(7,1) * t224 - mrSges(7,2) * t207) * t153 + (-mrSges(7,1) * t225 - mrSges(7,2) * t205) * t148)) * pkin(2); (-t186 - qJD(6) * (mrSges(7,1) * t253 + Ifges(7,5))) * t152 - t147 * (mrSges(7,1) * t196 - qJD(6) * (mrSges(7,2) * t253 + Ifges(7,6))); -qJD(6) * t181; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11), t13(16); t13(2), t13(3), t13(5), t13(8), t13(12), t13(17); t13(4), t13(5), t13(6), t13(9), t13(13), t13(18); t13(7), t13(8), t13(9), t13(10), t13(14), t13(19); t13(11), t13(12), t13(13), t13(14), t13(15), t13(20); t13(16), t13(17), t13(18), t13(19), t13(20), t13(21);];
Mq = res;
